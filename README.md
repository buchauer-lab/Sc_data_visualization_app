## Project introduction
This project is an interactive web application for visualizing and analyzing single cells data and their genes within different cell groups. It has been built in Python using the Shiny framework and relies on scanpy and the Anndata object. 
The app was developed  in the theoretical research group working in computational immunology leads by Lisa Buchauer, Professor of Systems Biology of Infectious Diseases at the Department of Infectious Diseases and Intensive Care at Charité - University Medicine, Berlin. 
See below for a installation instructions.


## Installation
First you need a working installation of Python 3.6 or later. If not, you can install Miniconda following the instructions of this website (https://docs.anaconda.com/miniconda/ ).
Then, to access the dataset, you will need to install Scanpy : 
conda install -c conda-forge scanpy python-igraph leidenalg
and pull scanpy from PyPI :
pip install scanpy
If you want to use PyPI only or having trouble installing Miniconda on Linux or Mac, please refer to this website : https://scanpy.readthedocs.io/en/stable/installation.html .

Next step is to install shiny :
pip install shiny


## Dataset Prepapration
To use the dataset, you need to do a preprocessing step on your data, so it comes as an anndata object with x,y,z precalculated features.

If you are using scanpy with anndata, please follow the following tutorials on scanpy :
-	Preprocessing and clustering 3k PBMCs (legacy workflow) : https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html 
-	Core plotting functions : https://scanpy.readthedocs.io/en/stable/tutorials/plotting/core.html  


## Interaction between different files
The file named “app” is the source code of the app. The file named “preprocessing_clustering” is the preprocessing and clustering tutorial. This file is useful because it provides the adata object the necessary information for the application.

The file named “core_plotting_function” is the core plotting function tutorial.



## Contribution :
Code contribution is welcomed especially if you know how to code the download buttons. I couldn’t figure out how to do it. Shiny documentation is available about download buttons (https://shiny.posit.co/py/api/core/ui.download_button.html#shiny.ui.download_button ).

Purpose of Download buttons : 
"Download table of downregulated genes" : " table in CSV format should be downloadable when the button is clicked. It contains the name, fold change value, and p-value of all the genes that appear in blue on the chart.

"Download table of upregulated genes" : " table in CSV format should be downloadable when the button is clicked. It contains the name, fold change value, and p-value of all the genes that appear in red on the chart.

## License

Scanpy :
SCANPY: large-scale single-cell gene expression data analysis
F. Alexander Wolf, Philipp Angerer, Fabian J. Theis
Genome Biology 2018 Feb 06. doi: 10.1186/s13059-017-1382-0.

Anndata :
anndata: Annotated data
Isaac Virshup, Sergei Rybakov, Fabian J. Theis, Philipp Angerer, F. Alexander Wolf
bioRxiv 2021 Dec 19. doi: 10.1101/2021.12.16.473007.