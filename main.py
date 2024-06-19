import os
import scanpy as sc
import anndata as ad
import logging
import numpy as np
import pandas as pd

sc.settings.set_figure_params(dpi=100, facecolor="white")


def setup_logging(log_file):
    """
    Sets up logging configuration.

    :param: log_file: Path to the log file.
    """
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def validate_path(path: str) -> bool:
    return os.path.isdir(path)


def mtx_to_h5ad(input_path: str, output_path: str) -> str:
    """
    Converts a single 10x Genomics data folder to AnnData format and saves it as h5ad file.
    Annotates mitochondrial data and performs a quality check

    :param: input_path: Path to the folder containing the 10x Genomics data.
    :param: output_path: Path to the folder where the h5ad file and QC plot will be saved.
    :return: Path to the saved h5ad file.
    """

    if not validate_path(input_path):
        logging.error(f"Input path is not a valid directory: {input_path}")
        raise

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    try:
        file_name = os.path.basename(input_path.rstrip('/'))
        logging.info(f"Pre-processing file: {file_name}")
        print(f"Pre-processing file: {file_name}")

        # Annotate mitochondrial genes, perform quality control
        adata = sc.read_10x_mtx(input_path, var_names='gene_symbols', cache=True)
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)
        output_path = os.path.join(output_path, f"qc.h5ad")
        sc.write(output_path, adata)

        # Plot QC metrics
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True
        )
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

        logging.info(f"Successfully processed file: {file_name}")
        print(f"Successfully processed file: {file_name}")

        return output_path

    except Exception as e:
        logging.error(f"An error occurred while pre-processing {input_path}: {e}")


def get_ann(input_path: str) -> ad.AnnData:
    """
    Fetches AnnData from a specified folder

    :param: input_path
    :return: AnnData
    """
    try:
        logging.info(f"Retrieving anndata from: {input_path}")
        adata = ad.read_h5ad(input_path)
        logging.info(f"Successfully retrieved anndata from: {input_path}")
        print(f"Retrieving anndata from: {input_path}")
        return adata
    except OSError as e:
        logging.error(f"An error occurred while retrieving {input_path}: {e}")
        raise
    except Exception as e:
        logging.error(f"An unexpected error occurred while retrieving {input_path}: {e}")
        raise


def check_ann(adata: ad.AnnData):
    """
    Gives basic information about anndata file

    :param: adata
    """
    try:
        logging.info(f"Retrieving information about .h5ad file")

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
        logging.error(f"An error occurred when retrieving information from .h5ad file: {e}")


def filter(adata, output_path) -> str:
    """
    Filters the AnnData object based on given criteria and saves the filtered object.

    :param: adata: AnnData object.
    :param: output_path: Path to save the filtered AnnData object.
    :return: Path to the saved filtered AnnData object.
    """
    min_cell = int(input("Minimum number of cells expressed: "))
    max_cell = int(input("Maximum number of cells expressed: "))
    min_genes = int(input("Minimum of genes per cell: "))
    max_genes = int(input("Maximum number of genes per cell: "))
    mito = int(input("Count mito: "))

    try:
        logging.info(f"Filtering AnnData")
        print(f"Filtering AnnData")

        output_path = os.path.join(output_path, f"filtered.h5ad")

        if min_cell > 0:
            sc.pp.filter_genes(adata, min_cells=min_cell)
            logging.info(f"Filtered genes that appear in less then {min_cell}")

        if max_cell > 0:
            sc.pp.filter_genes(adata, max_cells=max_cell)
            logging.info(f"Filtered genes that appear in more then {max_cell}")

        if min_genes > 0:
            sc.pp.filter_cells(adata, min_genes=min_genes)
            logging.info(f"Filtered cells with less then {min_genes}")

        if max_genes > 0:
            sc.pp.filter_cells(adata, max_genes=max_genes)
            logging.info(f"Filtered cells with more then {max_genes}")

        if mito > 0:
            if 'pct_counts_mt' not in adata.obs.columns:
                raise ValueError("'pct_counts_mt' column is missing in adata.obs")

            adata = adata[adata.obs['pct_counts_mt'] <= mito].copy()
            logging.info(f"Filtered cells with less then {mito}% mitchonrial genes")

        if os.path.exists(output_path):
            base_name, ext = os.path.splitext(output_path)
            count = 1
            while os.path.exists(f"{base_name}_{count}{ext}"):
                count += 1
            output_path = f"{base_name}_{count}{ext}"
            sc.write(output_path, adata)
        else:
            sc.write(output_path, adata)

        logging.info(f"AnnData filtered successfully")
        print(f"AnnData filtered successfully")
        return output_path

    except Exception as e:
        logging.error(f"An error occurred while filtering AnnData: {e}")
        raise


def normalize(adata, output_path):
    """
    Normalizes, Logarithmizes detects doublets and highly variable genes the AnnData object and saves it.

    :param adata: AnnData object.
    :param output_path: Path to save the normalized AnnData object.
    :return: Path to the saved normalized AnnData object.
    """
    try:
        logging.info(f"Normalizing AnnData")
        print(f"Normalizing AnnData")
        # Detecting doublets
        sc.pp.scrublet(adata)
        # Saving count data
        adata.layers["counts"] = adata.X.copy()
        # Normalizing to median total counts
        sc.pp.normalize_total(adata)
        # Logarithmize the data
        sc.pp.log1p(adata)
        # Marking highly variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.pl.highly_variable_genes(adata)

        output_path = os.path.join(output_path, f"normalized.h5ad")
        sc.write(output_path, adata)

        logging.info(f"AnnData normalized successfully")
        print(f"AnnData normalized successfully")
        return output_path

    except Exception as e:
        logging.error(f"An error occurred while normalizing AnnData: {e}")
        raise


def reduce_dimensions(adata, output_path):
    """
    Reduces dimensions of the AnnData object and saves it.

    :param: adata: AnnData object.
    :param: output_path: Path to save the reduced dimension AnnData object.
    :return: Path to the saved reduced dimension AnnData object.
    """
    try:
        logging.info(f"Reducing dimensions of AnnData")
        print(f"Reducing dimensions of AnnData")

        sc.tl.pca(adata)
        sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.pl.umap(adata, size=2,)

        output_path = os.path.join(output_path, f"pca_umap.h5ad")
        sc.write(output_path, adata)

        logging.info(f"Dimensions reduced successfully")
        print(f"Dimensions reduced successfully")
        return output_path

    except Exception as e:
        logging.error(f"Reducing dimensions failed: {e}")
        raise


def clustering(adata, output_path, n_iterations=2, random_seed=42, resolution=1.0, method='louvain'):
    """
    Clusters the AnnData object and saves it.

    :param: adata: AnnData object.
    :param: output_path: Path to save the clustered AnnData object.
    :param: n_iterations: Number of iterations for the clustering algorithm.
    :param: random_seed: Random seed for reproducibility.
    :param: resolution: Resolution parameter for the clustering algorithm.
    :param: method: Clustering method ('leiden' or 'louvain').
    :return: Path to the saved clustered AnnData object.
    """
    try:
        logging.info("Attempting to cluster adata")
        print("Attempting to cluster adata")

        np.random.seed(random_seed)
        sc.settings.seed = random_seed

        # Perform clustering
        if method == 'leiden':
            logging.info("Starting Leiden clustering.")
            # noinspection PyTypeChecker
            sc.tl.leiden(adata, flavor="igraph", n_iterations=n_iterations,
                         random_state=random_seed, resolution=resolution)
            logging.info(f"Leiden clustering completed with {n_iterations} iterations and resolution {resolution}")

        elif method == 'louvain':
            logging.info("Starting Louvain clustering.")
            sc.tl.louvain(adata, flavor="vtraag", random_state=random_seed, resolution=resolution)
            logging.info(f"Louvain clustering completed with resolution {resolution}")
        else:
            raise ValueError("Invalid method specified. Use 'leiden' or 'louvain'.")

        adata.uns['clustering_parameters'] = {
            'n_iterations': n_iterations,
            'random_seed': random_seed,
            'resolution': resolution,
            'clustering_method': method
        }

        # Plot UMAPs
        logging.info("Generating UMAP plots.")
        sc.pl.umap(adata, color=[method], legend_loc="on data")
        sc.pl.umap(adata, color=[method, "predicted_doublet", "doublet_score"], wspace=0.5, size=3,)
        sc.pl.umap(adata, color=[method, "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
                   wspace=0.5, ncols=2,)

        output_path = os.path.join(output_path, f"clustered_{method}.h5ad")

        if os.path.exists(output_path):
            base_name, ext = os.path.splitext(output_path)
            count = 1
            while os.path.exists(f"{base_name}_{count}{ext}"):
                count += 1
            output_path = f"{base_name}_{count}{ext}"

        sc.write(output_path, adata)
        logging.info(f"Clustering succeeded")
        print(f"Clustering succeeded")
        return output_path

    except Exception as e:
        logging.error(f"Clustering failed: {e}")
        raise


def mark_genes(adata, output_path):
    """
    Marks specific genes in the AnnData object.

    :param: adata AnnData object.
    :param: output_path Path to save outputs.
    """
    try:
        logging.info("Marking genes in AnnData")
        print("Marking genes in AnnData")

        # Perform ranking of genes based on groups
        sc.tl.rank_genes_groups(adata, groupby="louvain", method="wilcoxon")

        # Display a dot plot of ranked genes
        sc.pl.rank_genes_groups_dotplot(adata, groupby="louvain", standard_scale="var", n_genes=5)

        # Display UMAP plot colored by specific genes and louvain clusters
        sc.pl.umap(
            adata,
            color=["louvain", "SLPI", "ELANE", "COL3A1", "FBN1", "LUM"],
            frameon=False,
            ncols=3,
        )

        # Prompt for group number for further analysis
        group = input("Enter group number for further analysis (type END to end the program): ")

        while group != "END":
            # Get top 5 ranked genes for the specified group
            ranked_genes_df = sc.get.rank_genes_groups_df(adata, group=group)
            dc_cluster_genes = ranked_genes_df.head(5)["names"].tolist()

            # Update UMAP plot with top ranked genes for the specified group
            sc.pl.umap(
                adata,
                color=["louvain"] + dc_cluster_genes,  # Include louvain and top genes for color
                frameon=False,
                ncols=3,
            )

            # Prompt again for group number
            group = input("Enter group number for further analysis (type END to end the program): ")

        logging.info("Genes marked successfully")
        print("Genes marked successfully")

        try:
            logging.info("Saving ranked genes to file")
            print("Saving ranked genes to file")

            result = adata.uns['rank_genes_groups']
            groups = result['names'].dtype.names
            df = pd.DataFrame({group + '_' + key[:1]: result[key][group]
                               for group in groups for key in ['names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']})

            output_path = os.path.join(output_path, f"ranked_genes.csv")

            if os.path.exists(output_path):
                base_name, ext = os.path.splitext(output_path)
                count = 1
                while os.path.exists(f"{base_name}_{count}{ext}"):
                    count += 1
                output_path = f"{base_name}_{count}{ext}"

            df.to_csv(output_path)

            print(f"Ranked genes exported to {output_path}")
            logging.info("Ranked genes exported")

        except Exception as e:
            logging.error(f"Saving ranked genes to excel file failed: {e}")

    except Exception as e:
        logging.error(f"Marking failed: {e}")


def run_analysis(input_path, output_path):
    current_file = mtx_to_h5ad(input_path, output_path)
    adata = get_ann(current_file)
    current_file = filter(adata, output_path)
    adata = get_ann(current_file)
    current_file = normalize(adata, output_path)
    adata = get_ann(current_file)
    current_file = reduce_dimensions(adata, output_path)
    adata = get_ann(current_file)
    current_file = clustering(adata, output_path)
    adata = get_ann(current_file)
    mark_genes(adata, output_path)


if __name__ == "__main__":
    setup_logging("C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\mainlog.log")
    run_analysis("C:\\Users\\User\\Desktop\\pythonProject1\\testcase\\test1", "C:\\Users\\User\\Desktop\\pythonProject1\\rescase\\test1")
