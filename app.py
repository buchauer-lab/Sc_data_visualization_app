import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from shiny import App, render, ui, reactive
import scanpy as sc
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from shinywidgets import output_widget, render_widget
import io 
import csv

#get the data needed
adata = sc.read('.\write\pbmc3k_withoutX.h5ad')
adataX = sc.read('.\write\pbmc3k.h5ad')
pbmc = sc.datasets.pbmc68k_reduced()

marker_genes = [
    *["IL7R", "CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "CD14"],
    *["LGALS3", "S100A8", "GNLY", "NKG7", "KLRB1"],
    *["FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP"],
]

#----------------------------------------------------------
#code needed to plot the VOLCANO PLOTS

# Assuming cell types are annotated in 'cell_type' column in adata.obs / Perform differential expression analysis
sc.tl.rank_genes_groups(adataX, 'cell type', method='wilcoxon')

# Extract results
result = adataX.uns['rank_genes_groups']

# Initialize DataFrame to store results
results_list = []
groups = result['names'].dtype.names

for group in groups:
    # Get gene names, p-values, and adjusted p-values
    genes = result['names'][group]
    pvals = result['pvals'][group]
    pvals_adj = result['pvals_adj'][group]
    fold_changes = result['logfoldchanges'][group]
    ratio = np.exp2(fold_changes)
    negative_log_pval = -(np.log10(pvals_adj))

    # Create a DataFrame for this group
    group_df = pd.DataFrame({
        'gene': genes,
        'pval': pvals,
        'pval_adj': pvals_adj,
        'fold_change': fold_changes,
        'ratio' : ratio,
        'cell type': group,
        'negative_log_pval' : negative_log_pval
    })
    results_list.append(group_df)

# Concatenate all group results into a single DataFrame
results_df = pd.concat(results_list, ignore_index=True)

# Pivot the DataFrame to create separate columns for each cell type
results_pivot = results_df.pivot(columns='cell type', values=['gene', 'pval', 'pval_adj', 'fold_change', 'ratio', 'negative_log_pval'])

# Rename the columns
results_pivot.columns = [f'{col}_{cell_type}' for col, cell_type in results_pivot.columns]

# creation of upregulated genes dataframe
def create_upregulated_genes_df(cell_type, fold_change_threshold, pval_adj_threshold):
    fold_change_col = f'fold_change_{cell_type}'
    negative_log_pval_col = f'negative_log_pval_{cell_type}'
    gene_col = f'gene_{cell_type}'
    
    filtered_df = results_pivot.dropna(subset=[fold_change_col, negative_log_pval_col, gene_col])
    upregulated_genes = filtered_df[
        (filtered_df[fold_change_col] > fold_change_threshold) & 
        (filtered_df[negative_log_pval_col] > -np.log10(pval_adj_threshold))
    ]
    
    upregulated_genes_df = upregulated_genes[[gene_col, fold_change_col, negative_log_pval_col]].copy()
    upregulated_genes_df.columns = ['Gene', 'Fold Change', 'P-value (-log10)']
    
    return upregulated_genes_df

# creation of downregulated genes dataframe
def create_downregulated_genes_df(cell_type, fold_change_threshold, pval_adj_threshold):
    fold_change_col = f'fold_change_{cell_type}'
    negative_log_pval_col = f'negative_log_pval_{cell_type}'
    gene_col = f'gene_{cell_type}'
    
    filtered_df = results_pivot.dropna(subset=[fold_change_col, negative_log_pval_col, gene_col])
    downregulated_genes = filtered_df[
        (filtered_df[fold_change_col] < -fold_change_threshold) & 
        (filtered_df[negative_log_pval_col] > -np.log10(pval_adj_threshold))
    ]
    
    downregulated_genes_df = downregulated_genes[[gene_col, fold_change_col, negative_log_pval_col]].copy()
    downregulated_genes_df.columns = ['Gene', 'Fold Change', 'P-value (-log10)']
    
    return downregulated_genes_df
#-------------------------------------------------------------------
#Code part 1 : define the app UI
app_ui = ui.page_fluid(
    
    #Layout Part 1 : selection buttons
    ui.card(    #card 1 : contains the 3 cards for option selection     
        ui.card_header("UMAP AND VIOLIN PLOT"),
        
            ui.layout_columns(  #separation in 3 cards for options selection
            
            ui.card(    #card 1.1 option selection UMAP 1
                ui.input_selectize(  
                    "choiceumap",  
                    "Choose annotation for UMAP1:",  
                    {   "Cell type Leiden": "Cell type Leiden", 
                        "Cell type Louvain": "Cell type Louvain",
                        },  
                )  
            ),
            
            ui.card(    #card 1.2 option selection UMAP 2 and violin plot
                ui.input_selectize(
                    "choiceviolin",
                    "Choose value for UMAP2 and violin plot:",
                    {
                        "Gene expression": {"Gene expression": "Gene expression"},
                        " QC metrics": {
                            # "% mt reads": "% mt reads", 
                            "percent mito":"percent mito", 
                            "n genes" : "n genes"},
                    },
                ),
            ),
            
            ui.card(    #card 1.3 : Text field for Gene Name
                ui.input_text("textgene", "Gene name", ""),
                ui.output_text_verbatim("value_gene"),
            )
        ),
    ),
  
    
    #Layout Part 2 : UMAP + Violinplot
    ui.layout_columns(  #Separation between the 2 UMAP (on the left) and the Violin plot (on the right)
        
        ui.card(    #card 2.1 : contains 2 cards, one for each UMAP
            ui.card_header("card1 : UMAP"),
                                             
                ui.card(    #card umap 1 
                    ui.card_header("Comparative UMAP"),
                    ui.output_plot("umap_plot_comp", height = "350px", width="450px"),
                ),
             
            ui.card(    #card umap 2
                ui.card_header("UMAP of gene expression"),
                ui.output_plot("umap_plot_gene", height = "350px", width="450px"),
            )
        ),
        
        ui.card(    #card 2.2 violinplot
            ui.card_header("card 2 : Violin plots"),
            ui.card(
                ui.output_plot("violin_plot"),
                full_screen= True
            )                             
        )
    ),
    
    #Layout Part 3 : Volcano plot
    ui.card(
        ui.layout_sidebar(
            
            ui.sidebar(     #Creation of a sidebar, contains selection of cell type and text field to choose threeshold
                ui.input_selectize(     #selection
                        "choicevolcano",  
                        "Choose cell category",  
                        {"T-cell": "T-cell", "B-cell": "B-cell", "Monocytes": "Monocytes", "NK": "NK", "Monocytes": "Monocytes", "Dendritic": "Dendritic", "Plasma": "Plasma"},  
                    ),
                
                ui.input_text("textpvalue", "P-value threshold", "0.05"),   #textfield for p-value threshold
                ui.output_text_verbatim("value_pvalue"),
                
                ui.input_text("textfoldchange", "Fold change threshold", "1.5"),   #textfield for fold change threshold
                ui.output_text_verbatim("text_foldchange")
            ),

        output_widget("volcanoplot"),   #volcanoplot
        
        ui.layout_columns(      #separation in 2 part below the volcanoplot to add 2 buttons
            ui.card(
                ui.download_button("Download_downregulated","Download table of downregulated genes", class_="btn-primary")  #class_ : blue color
            ),
            ui.card(
                ui.download_button("Download_upregulated","Download table of upregulated genes", class_="btn-danger")   #class_ : red color
            )
        )        
        )
)
)

#---------------------------------------------------------------------------
#Code Part 2 : server logic
def server(input, output, session):
   
   #PART 1 : Getting in real time the user's selections
   @reactive.calc   #selection for UMAP1 and violin plot
   def choice_umap():
       return input.choiceumap
   
   @reactive.calc   #selection for UMAP2 and violin plot
   def choice_umap2_violin():
       return input.choiceviolin
   
   @reactive.calc   #text for gene name, umap 2 and violin plot
   def gene_exp():
        return input.textgene()
    
   @reactive.calc   #cell type selection for volcanoplot
   def choice_volcano():
       return input.choicevolcano
   
   @reactive.calc   #foldchange threshold value, volcanoplot
   def valuefoldchange():
       return float(input.textfoldchange())
   
   @reactive.calc   #p-value threshold value, volcanoplot
   def p_value():
       return float(input.textpvalue())
   

   #PART 2 : PLOTTING UMAP 1
   @render.plot
   def umap_plot_comp():
       choice_umap1 = choice_umap().get() #get selection for UMAP 1   

       if choice_umap1 == "Cell type Leiden" :
            ax = sc.pl.umap(adata, color="leiden", legend_loc="on data", title="UMAP cell type", frameon=False, save=".pdf", show=False)
            return ax
                
       elif choice_umap1 == "Cell type Louvain" :
            ax = sc.pl.umap(adataX,color="cell type",legend_loc="on data", frameon=False,legend_fontsize=10,legend_fontoutline=2, s = 75,show=False) 
            return ax
    
   #PART 3 : PLOTTING UMAP 2
   @render.plot
   def umap_plot_gene():
       choice_umap2 = choice_umap2_violin().get() #get selection for UMAP 2
       gene_name = gene_exp()   #get gene name 
       
       if choice_umap2 == "Gene expression" and gene_name :    #UMAP2 of gene expression 
            with rc_context({"figure.figsize": (3, 3)}):
                ax = sc.pl.umap(adata, color=gene_name, s=50, frameon=False, ncols=4, vmax="p99", show=False)
       
       elif choice_umap2 == "percent mito":    #UMAP2 of percent mito
            with rc_context({"figure.figsize": (3, 3)}):
                ax = sc.pl.umap(adataX, color="percent_mito", s=75, frameon=False, ncols=4, vmax="p99", show=False)
            return ax 
        
       elif choice_umap2 == "n genes":    #UMAP2 of n genes
            with rc_context({"figure.figsize": (3, 3)}):
                ax = sc.pl.umap(adataX, color= "n_counts", s=75, frameon=False, ncols=4, vmax="p99", show=False)
            
   #PART 4 : PLOTTING VIOLINPLOT
   @render.plot    
   def violin_plot():
        choice_x_axis = choice_umap().get()             #get selection for x axis
        choice_y_axis = choice_umap2_violin().get()     #get selection for y axis
        gene_name = gene_exp()                          #get gene name
           
        if choice_y_axis == "Gene expression" and  gene_name and choice_x_axis=="Cell type Leiden" :
            fig, ax = plt.subplots()
            sc.pl.violin(adata, gene_name, groupby="leiden", ax=ax, show=False)
            fig.set_size_inches(20, 3)  # Set the size of the figure
            
            # Rotate x-axis labels
            for label in ax.get_xticklabels():
                label.set_rotation(45)
                label.set_ha('right')  # Aligns the labels to the right side    
            return fig
        
        elif choice_y_axis == "Gene expression" and  gene_name and choice_x_axis=="Cell type Louvain":
            fig, ax = plt.subplots()
            sc.pl.violin(adataX, gene_name, groupby="cell type", ax=ax, show=False)
            fig.set_size_inches(20, 3)  # Set the size of the figure
            
            # Rotate x-axis labels
            for label in ax.get_xticklabels():
                label.set_rotation(45)
                label.set_ha('right')  # Aligns the labels to the right side   
            return fig
        
        elif choice_y_axis == "n genes" and choice_x_axis=="Cell type Leiden":
            fig, ax = plt.subplots()
            sc.pl.violin(
                adata,
                "n_genes",
                groupby="leiden",
                stripplot=False,  # remove the internal dots
                inner="box",  # adds a boxplot inside violins
                ax=ax,
                show=False
            )
            
            # Rotate x-axis labels
            for label in ax.get_xticklabels():
                label.set_rotation(45)
                label.set_ha('right')  # Aligns the labels to the right side
            return fig
        
        elif choice_y_axis == "n genes" and choice_x_axis=="Cell type Louvain":
            fig, ax = plt.subplots()
            sc.pl.violin(
                adataX,
                "n_genes",
                groupby="cell type",
                stripplot=False,  # remove the internal dots
                inner="box",  # adds a boxplot inside violins
                ax=ax,
                show=False
            )
            # Rotate x-axis labels
            for label in ax.get_xticklabels():
                label.set_rotation(45)
                label.set_ha('right')  # Aligns the labels to the right side
            return fig
        
        elif choice_y_axis == "percent mito" and choice_x_axis=="Cell type Leiden":
            fig, ax = plt.subplots()
            sc.pl.violin(
                adata,
                "pct_counts_mt",
                groupby="leiden",
                stripplot=False,  # remove the internal dots
                inner="box",  # adds a boxplot inside violins
                ax=ax,
                show=False
            )
            
            # Rotate x-axis labels
            for label in ax.get_xticklabels():
                label.set_rotation(45)
                label.set_ha('right')  # Optional: aligns the labels to the right side
                
            return fig
        
        elif choice_y_axis == "percent mito" and choice_x_axis=="Cell type Louvain":
            fig, ax = plt.subplots()
            sc.pl.violin(
                adataX,
                "percent_mito",
                groupby="cell type",
                stripplot=False,  # remove the internal dots
                inner="box",  # adds a boxplot inside violins
                ax=ax,
                show=False
            )
            
            # Rotate x-axis labels
            for label in ax.get_xticklabels():
                label.set_rotation(45)
                label.set_ha('right')  # Optional: aligns the labels to the right side
                
            return fig
                    
            
   #PART 4 : PLOTTING VOLCANOPLOT         
   @render_widget
   def volcanoplot():
        cell_type_volcano = choice_volcano().get()  #get cell type selection for volcanoplot

        fold_change_threshold = valuefoldchange()  #get foldchange threshold value
        pval_adj_threshold = p_value()             #get p-value threshold value
        
        #list of cell types, useful for the for loop
        cell_types = ["B-cell", "T-cell", "Monocytes","NK", "Dendritic", "Plasma"] 
        
        #calculation of values and posting of volcanoplot
        for cell_type in cell_types:
            
            if cell_type_volcano == cell_type:
                
                fold_change_col = f'fold_change_{cell_type}'
                negative_log_pval_col = f'negative_log_pval_{cell_type}'
                gene_col = f'gene_{cell_type}'
                
                fold_change = pd.to_numeric(results_pivot[fold_change_col], errors='coerce')
                negative_log_pval = pd.to_numeric(results_pivot[negative_log_pval_col], errors='coerce')
                                
                # Replace problematic values with NaN
                fold_change.replace([np.inf, -np.inf], np.nan, inplace=True)
                negative_log_pval.replace([np.inf, -np.inf], np.nan, inplace=True)
                    
                # Filter out rows with NaN values
                results_pivot[fold_change_col] = fold_change
                results_pivot[negative_log_pval_col] = negative_log_pval
                filtered_df = results_pivot.dropna(subset=[fold_change_col, negative_log_pval_col, gene_col])
            
                # Define the color based on the thresholds
                filtered_df['color'] = np.where(
                    (filtered_df[fold_change_col] > fold_change_threshold) & 
                    (filtered_df[negative_log_pval_col] > -np.log10(pval_adj_threshold)),
                    'Upregulated', np.where((filtered_df[fold_change_col] < -fold_change_threshold) & 
                            (filtered_df[negative_log_pval_col] > -np.log10(pval_adj_threshold)),
                            'Downregulated', 'Not significant') )
                
                # Identify the 5 most significant upregulated genes to write their names
                significant_up = filtered_df[filtered_df['color'] == 'Upregulated'].nlargest(5, negative_log_pval_col)
                
                #plotting the volcanoplot
                fig = px.scatter(data_frame=filtered_df, x=fold_change_col, y=negative_log_pval_col,
                        color="color", color_discrete_map={'Upregulated': 'red', 'Downregulated': 'blue',
                                            'Not significant': 'gray'},
                        title=f'Volcano Plot for {cell_type}', 
                        labels={fold_change_col: 'Fold Change', negative_log_pval_col: '-log10(p-value)'},
                        hover_data={gene_col: True, fold_change_col: True,
                                    negative_log_pval_col: True, 'color': False})
                
                # Update layout to customize the legend title
                fig.update_layout(legend_title_text='Regulation Status')
                
                # Update hover template to show gene_name first
                fig.update_traces(
                    hovertemplate=f'<b>%{{customdata[0]}}</b><br>Fold Change: %{{x}}<br>-log10(p-value): %{{y}}'
                    )
                
                # Add annotations for the most significant genes
                for i, row in significant_up.iterrows():
                    fig.add_annotation(x=row[fold_change_col], y=row[negative_log_pval_col],
                        text=row[gene_col], showarrow=True, arrowhead=1, ax=20, ay=-30)
                return fig
            
#run app     
app = App(app_ui, server, debug=True)