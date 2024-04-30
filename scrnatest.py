import os
import scanpy as sc
import anndata

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
    adata = anndata.read_h5ad(input_folder)
    return adata

def check_ann(adata):
    n_obs = adata.shape[0]
    n_var = adata.shape[1]
    print("Number of observations (cells):", n_obs)
    print("Number of variables (genes):", n_var)
    print("Available variables (annotations):", adata.var.keys())

def filter_ann(adata, output_folder,filter):
    sc.pp.filter_genes(adata, min_cells=filter)
    adata_filtered = adata[:, adata.obs['n_cells'] >= filter]
    output_path = output_folder + "filtered_data.h5ad"
    sc.write(output_path, adata_filtered)

i="C:\\Users\\User\\Desktop\\pythonProject1\\testcase"
o="C:\\Users\\User\\Desktop\\pythonProject1\\rescase"
ii="C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1.h5ad"

#mtx_to_h5ad(i,o)
adata=get_ann(ii)
check_ann(adata)
#filter_ann(adata,o,3)

