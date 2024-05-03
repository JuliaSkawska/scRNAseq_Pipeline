import os
import re
import scanpy as sc
import anndata
import pandas as pd
import logging

def setup_logging(log_file):
    '''
    Set up logging configuration.

    :param log_file: Path to the log file.
    '''
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
def mtx_to_h5ad(input_folder, output_folder):
    '''
    Convert 10x Genomics data to AnnData format and save as h5ad files.

    :param input_folder: Path to the folder containing 10x Genomics files.
    :param output_folder: Path to the folder where the h5ad files will be saved.
    '''
    for file_name in os.listdir(input_folder):
        folder_path = os.path.join(input_folder, file_name)
        if os.path.isdir(folder_path):
            try:
                logging.info(f"Processing folder from mxt_to_h5ad: {file_name}")

                adata = sc.read_10x_mtx(folder_path, var_names='gene_symbols', cache=True)

                output_folder_path = os.path.join(output_folder, file_name)
                os.makedirs(output_folder_path, exist_ok=True)
                output_file = os.path.join(output_folder_path, f"{file_name}.h5ad")
                adata.write_h5ad(output_file)

                logging.info(f"Successfully processed folder from mxt_to_h5ad: {file_name}")
            except Exception as e:
                logging.error(f"An error occurred while processing {file_name}: {e}")

def get_ann(input_folder):
    '''
    returns anndata from a specified folder
    :param input_folder:
    :return:
    '''
    try:
        adata = anndata.read_h5ad(input_folder)
        logging.info(f"Successfully retrived anndata: {input_folder}")
        return adata
    except Exception as e:
        print(f"An error occurred while reading the AnnData file: {e}")
        return None

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

if __name__ == "__main__":
    i = "C:\\Users\\User\\Desktop\\pythonProject1\\testcase"
    o = "C:\\Users\\User\\Desktop\\pythonProject1\\rescase"
    ii = "C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\filtered_data.h5ad"
    iii = "C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\annotated_ann.h5ad"
    iiii = "C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\qc.h5ad"

    setup_logging("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\mainlog.log")
    #mtx_to_h5ad(i, o)

    adata=get_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\test1.h5ad")
    #quality_check(adata,o)
    #adata=get_ann(iiii)
    #check_ann(adata)
    #annotate_ann(adata,o)
    #adata=get_ann(iii)
    #check_ann(adata)


