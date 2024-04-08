#import os
import scanpy as sc
#import pandas as pd


# Creating and saving anndata from matrix
def mtx_to_h5ad(input_path, output_path):
    adata = sc.read_10x_mtx(input_path,var_names='gene_symbols',cache=True)

    adata.write_h5ad(output_path)

input_path = input("Input the path to the directory containing the scRNA-seq data: ")
output_path = input("Input the path to the output file (including filename and extension): ")

mtx_to_h5ad(input_path, output_path)
