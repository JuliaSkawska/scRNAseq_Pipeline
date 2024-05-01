import os
import re
import scanpy as sc
import anndata
import pandas as pd

def mtx_to_h5ad(input_folder, output_folder):
    '''
    Convert 10x Genomics data to AnnData format and save as h5ad files.

    :param input_folder: Path to the folder containing 10x Genomics files.
    :param output_folder: Path to the folder where the h5ad files will be saved.
    '''
    for folder_name in os.listdir(input_folder):
        folder_path = os.path.join(input_folder, folder_name)
        if os.path.isdir(folder_path):
            adata = sc.read_10x_mtx(folder_path, var_names='gene_symbols', cache=True)
            output_filename = os.path.join(output_folder, f"{folder_name}.h5ad")
            adata.write_h5ad(output_filename)

def get_ann(input_folder):
    '''
    returns anndata from a specified folder
    :param input_folder:
    :return:
    '''
    adata = anndata.read_h5ad(input_folder)
    return adata

def check_ann(adata):
    '''
    Gives basic information about anndata file
    :param adata:
    :return:
    '''
    n_obs = adata.shape[0]
    n_var = adata.shape[1]
    print("Number of observations (cells):", n_obs)
    print("Number of variables (genes):", n_var)
    print("Available variables (annotations):", adata.var.keys())

def filter_ann(adata, output_folder,filter):
    '''
    Filters anndata file for minimal amount of cells specified in the filter param
    :param adata:
    :param output_folder:
    :param filter:
    :return:
    '''
    sc.pp.filter_genes(adata, min_cells=filter)
    output_path = output_folder + "\\filtered_data.h5ad"
    sc.write(output_path, adata)

def annotate_ann(adata, output_folder):
    '''
    Checks if gene in anndata file is mitochondrail, returns a text file with gene name and True/False statement
    :param adata:
    :param output_folder:
    '''
    gene_names = adata.var_names
    output_file = output_folder + "\\gene_annotation.txt"

    with open(output_file, 'w') as f:
        for gene_name in gene_names:
            if gene_name.startswith("MT"):
                f.write(gene_name + "\tTrue\n")
            else:
                f.write(gene_name + "\tFalse\n")


i="C:\\Users\\User\\Desktop\\pythonProject1\\testcase"
o="C:\\Users\\User\\Desktop\\pythonProject1\\rescase"
ii="C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\filtered_data.h5ad"

#mtx_to_h5ad(i,o)
adata=get_ann(ii)
#check_ann(adata)
#filter_ann(adata,o,3)
annotate_ann(adata,"C:\\Users\\User\\Desktop\\pythonProject1\\rescase")


