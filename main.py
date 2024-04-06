import scanpy as sc
import pandas as pd

#creating and saving anndata from matrix
adata = sc.read_10x_mtx('C:/Users/User/Desktop/pythonProject1/testcase/',
                      var_names='gene_symbols',
                      cache=True)

adata.write_h5ad("C:/Users/User/Desktop/pythonProject1/output_file.h5ad")



