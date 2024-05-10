import os
import re
import scanpy as sc
import anndata
import pandas as pd
import logging
import matplotlib.pyplot as plt

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
        logging.info(f"Retrieving anndata from: {input_folder}")
        adata = anndata.read_h5ad(input_folder)
        logging.info(f"Successfully retrieved anndata from: {input_folder}")
        return adata
    except Exception as e:
        logging.error(f"An error occurred while processing {input_folder}: {e}")
        return None

def check_ann(adata):
    '''
    Gives basic information about anndata file
    :param adata:
    :return:
    '''
    try:
        logging.info(f"Retrieving anndata info")
        n_obs = adata.shape[0]
        n_var = adata.shape[1]
        print("Number of observations (cells):", n_obs)
        print("Number of variables (genes):", n_var)
        if hasattr(adata, 'var') and hasattr(adata.var, 'keys'):
            print("Available variables (annotations):", adata.var.keys())
        else:
            print("No variable annotations available.")

        if hasattr(adata, 'obs'):
            print("Available observations (annotations):", adata.obs)
        else:
            print("No observation annotations available.")

    except Exception as e:
        logging.error(f"An error occurred when retrieving info from anndata: {e}")

def filter_ann(adata, output_folder,filter_value):

    try:
        if not isinstance(filter_value, int) or filter_value <= 0:
            raise ValueError("Filter value must be a positive integer.")

        if not os.path.isdir(output_folder):
            os.makedirs(output_folder)
            logging.info(f"Created output folder: {output_folder}")

        logging.info("Filtering Anndata...")
        sc.pp.filter_genes(adata, min_cells=filter_value)
        output_path = os.path.join(output_folder, "filtered_data.h5ad")
        sc.write(output_path, adata)
        logging.info("AnnData filtered successfully.")
    except Exception as e:
        logging.error(f"An error occurred while filtering AnnData: {e}")

def annotate_ann(adata,file_path):
    '''
    Checks if genes in anndata file are mitochondrial, adds a column named mito to AnnData object with True/False statement

    :param adata: AnnData object containing gene information.
    '''
    mito_column=[]
    try:
        logging.info("Attempting to annotate anndata with mito")
        for gene_name in adata.var_names:
            if re.match(r'[Mm][Tt]', gene_name):
                mito_column.append(True)
            else:
                mito_column.append(False)

        if len(mito_column) != len(adata.var):
            raise ValueError("Length of new_annotations does not match the number of genes in adata.var")

        adata.var['mito'] = mito_column
        adata.write(file_path+"\\annotated_ann.h5ad")
        logging.info("Anndata successfully annotated")
    except Exception as e:
        logging.error(f"An error occurred while annotating AnnData: {e}")

def filter_ann(adata, file_path, min_gene, max_gene, min_cell, max_cell, normalize=True, log_transform=True):
    '''
    Filters anndata file for minimal amount of cells specified in the filter param
    :param adata:
    :param file_path ( output folder ):
    :param min_gene:
    :param min_cell:
    :param normalize:
    :param log_transform:
    :return:
    '''
    try:
        logging.info("Attempting to filter AnnData")
        if min_gene>0:
            sc.pp.filter_cells(adata, min_genes=min_gene)
        if max_gene>0:
            sc.pp.filter_cells(adata, max_genes=max_gene)
        if min_cell>0:
            sc.pp.filter_genes(adata, min_cells=min_cell)
        if max_cell>0:
            sc.pp.filter_genes(adata, max_cells=max_cell)
        if normalize==True:
            sc.pp.normalize_total(adata)
        if log_transform==True:
            sc.pp.log1p(adata)
        adata.write(file_path + "\\filtered_ann.h5ad")
        logging.info("AnnData filtered successfully")
    except Exception as e:
        logging.error(f"An error occurred while filtering AnnData: {e}")

def quality_check(adata, file_path):
    try:
        logging.info("Attempting to calculate qc metrics for AnnData")
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        output_file = os.path.join(file_path, "qc.h5ad")
        adata.write(output_file)
        logging.info("Calculated qc metrics for AnnData and saved to file")
    except Exception as e:
        logging.error(f"An error occurred while calculating qc metrics for AnnData: {e}")

def violin_plot(adata,file_path):
    try:
        logging.info("Attempting to generate a plot for AnnData")
        sc.pl.violin(adata, keys=['n_genes_by_counts', 'total_counts', 'pct_counts_mito'], multi_panel=True)
        output_file = os.path.join(file_path, 'violin_plot.png')
        plt.savefig(output_file)
        logging.info("Violin plot was sucessfully generated")
    except Exception as e:
        logging.error(f"An error occured while trying to generate a AnnData plot: {e}")

if __name__ == "__main__":
    i = "C:\\Users\\User\\Desktop\\pythonProject1\\testcase"
    u = "C:\\Users\\User\\Desktop\\pythonProject1\\rescase"
    o = "C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1"
    ii = "C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\filtered_ann.h5ad"
    iii = "C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\annotated_ann.h5ad"
    iiii = "C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\qc.h5ad"
    setup_logging("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\mainlog.log")
    #mtx_to_h5ad -> filter -> annotate -> quality
    #mtx_to_h5ad(i, u)
    #adata=get_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\test1.h5ad")
    #check_ann(adata)
    #adata=get_ann(iii)
    #filter_ann(adata,o,200,6000,3,0,True,True)
    #check_ann(adata)
    #annotate_ann(adata,o)
    #adata=get_ann(ii)
    #check_ann(adata)
    #quality_check(adata,o)
    #check_ann(adata)
    #adata=get_ann(iiii)
    #violin_plot(adata,o)



