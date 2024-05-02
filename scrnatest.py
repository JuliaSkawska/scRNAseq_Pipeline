import os
import re
import scanpy as sc
import anndata
import pandas as pd
import logging

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
    print("Available observations ( annotations): ", adata.obs)

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

def annotate_ann(adata,file_path):
    '''
    Checks if genes in anndata file are mitochondrial, adds a column named mito to AnnData object with True/False statement

    :param adata: AnnData object containing gene information.
    '''
    mito_column=[]

    for gene_name in adata.var_names:
        if re.match(r'[Mm][Tt]', gene_name):
            mito_column.append(True)
        else:
            mito_column.append(False)

    if len(mito_column) != len(adata.var):
        raise ValueError("Length of new_annotations does not match the number of genes in adata.var")
    adata.var['mito'] = mito_column
    adata.write(file_path+"\\annotated_ann.h5ad")

i="C:\\Users\\User\\Desktop\\pythonProject1\\testcase"
o="C:\\Users\\User\\Desktop\\pythonProject1\\rescase"
ii="C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\filtered_data.h5ad"
iii="C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\annotated_ann.h5ad"
iiii="C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\qc.h5ad"


#adata=get_ann(iii)
#quality_check(adata,o)
#adata=get_ann(iiii)
#check_ann(adata)
#annotate_ann(adata,o)
#adata=get_ann(iii)
#check_ann(adata)


