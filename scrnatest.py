import os
import re
import scanpy as sc
import anndata
import pandas as pd
import logging
import numpy as np
import matplotlib.pyplot as plt


def setup_logging(log_file):
    """
    Sets up logging configuration.

    :param log_file: Path to the log file.
    """
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def validate_path(path):
    """
    Validates if the given path is a valid directory.

    :param path: Path to be validated.
    :return: True if the path is a valid directory, False otherwise.
    """
    return os.path.isdir(path)


def mtx_to_h5ad(input_path, output_path):
    """
    Converts 10x Genomics data to AnnData format and save as h5ad files.

    :param input_path:
    :param output_path: Path to the folder where the h5ad files will be saved
    """
    if not validate_path(input_path):
        logging.error(f"Input path is not a valid directory: {input_path}")
        return
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for file_name in os.listdir(input_path):
        folder_path = os.path.join(input_path, file_name)
        if os.path.isdir(folder_path):
            try:
                logging.info(f"Processing folder from mxt_to_h5ad: {file_name}")

                adata = sc.read_10x_mtx(folder_path, var_names='gene_symbols', cache=True)

                output_folder_path = os.path.join(output_path, file_name)
                os.makedirs(output_folder_path, exist_ok=True)
                output_file = os.path.join(output_folder_path, f"{file_name}.h5ad")
                adata.write_h5ad(output_file)

                logging.info(f"Successfully processed folder from mxt_to_h5ad: {file_name}")
            except Exception as e:
                logging.error(f"An error occurred while processing {file_name}: {e}")


def get_ann(input_path):
    """
    Fetches AnnData from a specified folder

    :param input_path:
    :return: AnnData
    """
    try:
        logging.info(f"Retrieving anndata from: {input_path}")
        adata = anndata.read_h5ad(input_path)
        logging.info(f"Successfully retrieved anndata from: {input_path}")
        return adata
    except OSError as e:
        logging.error(f"An error occurred while processing {input_path}: {e}")
        raise
    except Exception as e:
        logging.error(f"An unexpected error occurred while processing {input_path}: {e}")
        raise


def check_ann(adata):
    """
    Gives basic information about anndata file

    :param adata:
    """
    try:
        logging.info(f"Retrieving anndata information")

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


def filter_ann(output_path, adata, min_gene=None, max_gene=None, min_cell=None, max_cell=None, mito=None, normalize=True, log_transform=True):
    """
    Filters AnnData file for minimal amount of cells specified in the filter parameters.

    :param output_path: Path to save the output AnnData file
    :param adata: AnnData object
    :param min_gene: Minimum number of genes per cell
    :param max_gene: Maximum number of genes per cell
    :param min_cell: Minimum number of cells per gene
    :param max_cell: Maximum number of cells per gene
    :param normalize: Whether to perform total count normalization (default is True)
    :param log_transform: Whether to perform log transformation (default is True)
    """
    try:
        logging.info("Starting to filter AnnData")

        if min_gene is not None:
            sc.pp.filter_cells(adata, min_genes=min_gene)
            logging.info(f"Filtered cells with fewer than {min_gene} genes")
        if max_gene is not None:
            sc.pp.filter_cells(adata, max_genes=max_gene)
            logging.info(f"Filtered cells with more than {max_gene} genes")

        if min_cell is not None:
            sc.pp.filter_genes(adata, min_cells=min_cell)
            logging.info(f"Filtered genes with fewer than {min_cell} cells")
        if max_cell is not None:
            sc.pp.filter_genes(adata, max_cells=max_cell)
            logging.info(f"Filtered genes with more than {max_cell} cells")

        if normalize:
            sc.pp.normalize_total(adata)
            logging.info("Performed total count normalization")

        if log_transform:
            sc.pp.log1p(adata)
            logging.info("Performed log transformation")

        if mito is not None:
            if 'pct_counts_mito' not in adata.obs.columns:
                raise ValueError("'pct_counts_mito' column is missing in adata.obs")

            adata = adata[adata.obs['pct_counts_mito'] <= mito].copy()
            logging.info(f"Filtered observations with pct_counts_mito <= {mito}")


        filtered_file_path = os.path.join(output_path, "filtered_ann.h5ad")
        adata.write(filtered_file_path)
        logging.info(f"Filtered AnnData saved successfully to: {filtered_file_path}")

        return filtered_file_path

    except Exception as e:
        logging.error(f"Error occurred while filtering AnnData: {e}")
        raise


def annotate_ann(output_path, adata):
    """
    Checks if genes in anndata file are mitochondrial, adds a column named mito to AnnData object with bool value

    :param output_path:
    :param adata: AnnData object containing gene information.
    """
    mito_column = []
    try:
        logging.info("Attempting to annotate AnnData with mitochondrial genes")

        for gene_name in adata.var_names:
            if re.match(r'[Mm][Tt]', gene_name):
                mito_column.append(True)
            else:
                mito_column.append(False)

        if len(mito_column) != len(adata.var):
            raise ValueError("Length of mito_column does not match the number of genes in adata.var")

        adata.var['mito'] = mito_column

        annotated_file_path = os.path.join(output_path, "annotated_ann.h5ad")
        adata.write(annotated_file_path)

        logging.info(f"Annotated AnnData saved successfully to: {annotated_file_path}")
        return annotated_file_path

    except Exception as e:
        logging.error(f"An error occurred while annotating AnnData: {e}")


def quality_check(output_path, adata):
    """
    Calculates quality control (QC) metrics for AnnData.

    :param output_path: The path to save the QC metrics file.
    :param adata: The AnnData object containing the data.
    """
    try:
        logging.info("Calculating QC metrics for AnnData")

        sc.pp.calculate_qc_metrics(adata, inplace=True)
        qc_output_file = os.path.join(output_path, "qc.h5ad")
        adata.write(qc_output_file)

        logging.info("QC metrics calculated and saved to file: %s", qc_output_file)
        return qc_output_file

    except Exception as e:
        logging.error(f"Error calculating QC metrics for AnnData: {e}")

def add_pct_counts_mito(output_path, adata):
    """
    Adds pct_counts_mito metric to AnnData file (mitochondrial genes/all genes)
    :param output_path: Path to save the output AnnData file
    :param adata: AnnData object
    """
    try:
        logging.info("Attempting to add pct_counts_mito to AnnData file")

        if 'mito' not in adata.var.columns:
            raise ValueError("'mito' column is missing in adata.var")

        total_counts_per_cell = np.sum(adata.X, axis=1)

        # Boolean mask for mitochondrial genes
        mito_mask = adata.var['mito'].values

        mito_counts_per_cell = np.sum(adata.X[:, mito_mask], axis=1)

        pct_counts_mito = (mito_counts_per_cell / total_counts_per_cell) * 100

        adata.obs['pct_counts_mito'] = pct_counts_mito

        pct_file_path = os.path.join(output_path, "pct_mito.h5ad")
        adata.write(pct_file_path)
        logging.info("pct_counts_mito successfully added to file: %s", pct_file_path)

    except Exception as e:
        logging.error(f"Error when adding pct_counts_mito to AnnData: {e}")

def violin_plot(output_path, adata):
    """
    Generates a violin plot for AnnData.

    :param output_path: The path to save the violin plot image.
    :param adata: The AnnData object containing the data.
    """
    try:
        logging.info("Generating violin plot for AnnData")

        sc.pl.violin(adata, keys=['n_genes_by_counts', 'total_counts', 'pct_counts_mito'], multi_panel=True)
        plot_output_file = os.path.join(output_path, 'violin_plot.png')
        plt.savefig(plot_output_file)

        logging.info("Violin plot generated and saved to file: %s", plot_output_file)

    except Exception as e:
        logging.error(f"Error generating violin plot for AnnData: {e}")

if __name__ == "__main__":
    #mtx_to_h5ad -> filter (<3 cells ) -> annotate -> add pct_counts_mito  -> quality_check -> visualise -> filter -> visualise
    setup_logging("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\mainlog.log")
    op="C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1"
    #mtx_to_h5ad("C:\\Users\\User\\Desktop\\pythonProject1\\testcase","C:\\Users\\User\\Desktop\\pythonProject1\\rescase")
    #adata=get_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\test1.h5ad")
    #filter_ann(output_path=op,adata=adata,min_cell=3,log_transform=False,normalize=False)
    #adata=get_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\filtered_ann.h5ad")
    #annotate_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1",adata)
    #adata=get_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\annotated_ann.h5ad")
    #add_pct_counts_mito("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1",adata)
    #adata=get_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\pct_mito.h5ad")
    #quality_check("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1",adata)
    #adata=get_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\qc.h5ad")
    #violin_plot("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1",adata)
    #adata=get_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\qc.h5ad")
    #filter_ann(output_path=op,adata=adata,min_gene=200,max_gene=6000,mito=15,normalize=False,log_transform=False)
    #adata=get_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\filtered_ann.h5ad")
    #add_pct_counts_mito("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1", adata)
    #adata=get_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\pct_mito.h5ad")
    #quality_check("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1",adata)
    #adata = get_ann("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1\\qc.h5ad")
    #violin_plot("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1",adata)





