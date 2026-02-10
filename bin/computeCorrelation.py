import argparse
import scanpy as sc
import anndata as ad
import gc
from pathlib import Path
import numpy as np
import scipy.sparse
import pandas as pd 
from scipy.stats import spearmanr
import bottleneck 
from pathlib import Path
import bottleneck as bn
import psutil
import os 

def rank(data):
    """Rank normalize data with NaN handling"""
    orig_shape = data.shape
    data = bn.nanrankdata(data).reshape(-1) - 1  # Ranks, ignoring NaNs
    # Normalize between 0 and 1
    return (data / np.nansum(~np.isnan(data))).reshape(orig_shape)

def sparse_corr(A):
    """Compute Pearson correlation for a sparse matrix.
    
    Arguments:
        A {scipy.sparse matrix} -- Sparse matrix of gene expression values
    Returns:
        np.matrix -- Correlation matrix
    """
    dense_matrix = A.toarray()
    COR = np.corrcoef(dense_matrix, rowvar=False)
    return COR


def process_sample(sample_id, adata, min_cells_threshold, correlation_type, replace_nans):
    """
    Process data for a single sample: compute correlation and rank normalize.

    Returns:
        tuple -- Sample ID and a dictionary with processed edge lists and rank-normalized matrix.
    """
    subset_data = adata[adata.obs['sample_id'] == sample_id].X
    #subset_data = scipy.sparse.csr_matrix(subset_data)

    if subset_data.shape[0] < min_cells_threshold:
        print(f"Skipping sample {sample_id} due to insufficient number of cells.")
        return sample_id, None
    
    # Compute correlation matrix
    if correlation_type == 'pearson':
        corr_matrix = sparse_corr(scipy.sparse.csr_matrix(subset_data))
    elif correlation_type == 'spearman':
        corr_matrix, _ = spearmanr(subset_data, axis=0, nan_policy='omit')
    else:
        raise ValueError(f"Unsupported correlation type: {correlation_type}")
    
    np.fill_diagonal(corr_matrix, 1)
    rank_normalized_matrix = rank(corr_matrix) 

    if replace_nans:
        rank_normalized_matrix[np.isnan(rank_normalized_matrix)] = bottleneck.nanmean(rank_normalized_matrix)

    print(f"Finished processing sample {sample_id}")
    del corr_matrix 
    #del subset_data
    return sample_id, {
        'rank_normalized_matrix': rank_normalized_matrix,
        'gene_names': adata.var_names.tolist()
    }

def cor_mtx_to_edgeList(file_path, agg_adata, edge_list_dir2 ):

    # 1.  data into an edgeList
    df = pd.DataFrame(agg_adata.X, index=agg_adata.var_names, columns=agg_adata.var_names)
    edge_list = df.stack().reset_index()
    #edge_list.columns = ["source", "target", "weight"]
    edge_list.columns = ["geneA", "geneB", "rank_norm_corr"]
    
    # Remove self-loops & multiIndex
    edge_list_clean = edge_list[edge_list["geneA"] != edge_list["geneB"]]
    edge_list_clean.set_index(['geneA', 'geneB'], inplace=True)
    
    # Count the number of targets per source
    target_counts = edge_list_clean.groupby(level=0).size()
    target_counts_df = target_counts.reset_index()
    target_counts_df.columns = ["geneA", "num_targets"]
    
    # Sanity check
    if (target_counts_df["num_targets"] == (agg_adata.shape[0]-1)).all():
    #if (target_counts_df["num_targets"] == 10907).all():
        print(f"Sanity check passed")
        
        # Save edge list 
        output_file = edge_list_dir2 / f"{file_path.stem}_edgeList.csv"
        edge_list_clean.to_csv(output_file)
        print(f"Saved edge list to: {output_file}\n")
    else:
        print(f"Sanity check failed. Skipping file.\n")

def correlation_matrix_to_edgelist(corr_matrix: pd.DataFrame):

    upper_triangle_mask = np.triu(np.ones(corr_matrix.shape), k=1) ## mask the upper triangle excluding the diagonal
    
    # Apply the mask to the correlation matrix to get the non-redundant values
    upper_triangle_corr = corr_matrix.where(upper_triangle_mask == 1)
    
    # Unstack the upper triangle into a series, which will automatically drop NaN values
    edge_list = upper_triangle_corr.stack().reset_index()
    
    # Rename columns to match the desired edge list format
    edge_list.columns = ['geneA', 'geneB', 'rank_normalized_correlation']
    
    # Optionally, you can set geneA and geneB as the index
    edge_list.set_index(['geneA', 'geneB'], inplace=True)
    
    return edge_list

def save_to_hdf5(cor, output_dir, file_path):
    """
    Save correlation data, including edge lists and aggregated sparse matrices, to specified directories.
    
    Args:
        cor (dict): The correlation dictionary containing matrices, gene names, and data summary.
        output_dir (Path): The base output directory to save the data.
        file_path (Path): The input file path (used for naming output files).
    """
    # Create required directories
    agg_mtx_dir = output_dir / "agg-mtx"
    edge_list_dir1 = output_dir / "edge-list" / "upTriang"
    edge_list_dir2 = output_dir / "edge-list" / "all-matrix"
    summary_dir = output_dir / "cor-summary"

    agg_mtx_dir.mkdir(parents=True, exist_ok=True)
    edge_list_dir1.mkdir(parents=True, exist_ok=True)
    edge_list_dir2.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    # File names
    agg_file = agg_mtx_dir / f"{file_path.stem}_corAggSparse.h5ad"
    edge_list_file = edge_list_dir1 / f"{file_path.stem}_corEdgeLists.csv"
    #edge_list_file = edge_list_dir / f"{file_path.stem}_corEdgeLists.csv.gz"
    df_summary_file = summary_dir / f"{file_path.stem}_corSummary.csv"

    #1. Save AnnData object
    aggregated_matrix = cor["aggregated_rank_normalized_matrix"]

    gene_names = cor["gene_names"]
    agg_adata = ad.AnnData(X=aggregated_matrix)
    agg_adata.var_names = gene_names
    agg_adata.write_h5ad(agg_file, compression="gzip")
    print(f"Saved aggregated matrix and gene names to {agg_file}")

    #2. Save edgeList
    aggregated_matrixEdge = cor["edgeList"] 
    edge_list= correlation_matrix_to_edgelist(aggregated_matrixEdge)
    edge_list.to_csv(edge_list_file)
    print("saved edge list")

    #2.2 Save edgeList
    cor_mtx_to_edgeList(file_path, agg_adata, edge_list_dir2)

    #3. Save data summary 
    #cor["data_summary"].to_csv(df_summary_file, index=False)
    #print(f"Saved data summary to {df_summary_file}")

def compute_gene_correlation(adata, 
                             correlation_type='pearson', 
                             replace_nans=True, 
                             min_cells_threshold=20, 
                             file_path="", 
                             output_dir="",
                             xCellxSubject=False):
    """
    Compute gene-gene correlation matrices.
    
    Parameters
    ----------
    adata : AnnData
        Single-cell data.
    correlation_type : str
        'pearson' or 'spearman'.
    replace_nans : bool
        Whether to replace NaNs.
    min_cells_threshold : int
        Minimum number of cells to include a sample.
    category : str or None
        Disease or condition to subset. If 'ALL' or None, use all samples.
    xCellxSubject : bool
        If True, compute correlation across all samples at once; 
        if False, compute per-sample and aggregate.
    """
    print("Starting correlation computation...")
    file_path = Path(file_path)
    output_dir = Path(output_dir)

    subset = adata
    gene_names = subset.var_names
    aggregated_matrix = None
    sample_count = 0

    
    print(f"Computing correlation per sample.")
    samples_ids = subset.obs['sample_id'].unique()
    total_samples = len(samples_ids)
    print(f"  Total number of samples: {total_samples}")

    for sample in samples_ids:
        print(f"  Processing sample ID: {sample}")
        sample_id, result = process_sample(sample, subset, min_cells_threshold, correlation_type, replace_nans)

        if result is not None:
            sample_matrix = result["rank_normalized_matrix"]
            if aggregated_matrix is None:
                aggregated_matrix = sample_matrix.copy()
            else:
                aggregated_matrix += sample_matrix
            sample_count += 1
            print(f"  Sample {sample} matrix aggregated.")
            del sample_matrix
            gc.collect()
        else:
            print(f"  No valid matrix for sample {sample}.")

    if aggregated_matrix is not None and sample_count > 0:
        aggregated_matrix = rank(aggregated_matrix)

    # Convert aggregated matrix to DataFrame
    aggregated_matrixEdge = pd.DataFrame(aggregated_matrix, index=gene_names, columns=gene_names)

    # Summary
    summary_df = pd.DataFrame({
        "total_genes": [len(gene_names)],
        "total_single_cells": [subset.shape[0]],
        "total_samples": [sample_count]
    })

    category_data = {
        'aggregated_rank_normalized_matrix': aggregated_matrix,
        'edgeList': aggregated_matrixEdge, 
        'gene_names': gene_names,
        'data_summary': summary_df, 
        'file_path': str(file_path)
    }
    category_file_path = Path(category_data["file_path"])

    # Log memory usage
    process = psutil.Process(os.getpid())
    print(f" Memory usage before saving: {process.memory_info().rss / 1024 ** 2:.2f} MB")

    print(f" Saving results.")
    save_to_hdf5(category_data, output_dir, category_file_path)
    print(f" Results saved successfully to {category_file_path}.")

    # Clean up memory
    del subset, aggregated_matrix, aggregated_matrixEdge, category_data, summary_df
    gc.collect()
    print(f" Memory usage after cleanup: {process.memory_info().rss / 1024 ** 2:.2f} MB")


def main():
    parser = argparse.ArgumentParser(description="Compute gene correlation for a single file and category.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to input .h5ad file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save results.")
    parser.add_argument("--category", type=str, default=None, help="Disease category (or 'ALL' for all samples).")
    parser.add_argument("--correlation_type", type=str, default="pearson",
                        choices=["pearson", "spearman"], help="Type of correlation.")
    parser.add_argument("--replace_nans", action='store_true', help="Replace NaNs with median.")
    parser.add_argument("--xCellxSubject", action='store_true',
                        help="If set, compute correlation across all samples instead of per-sample aggregation.")
    parser.add_argument("--min_cells_threshold", type=int, default=20,
                        help="Minimum number of cells required per sample.")

    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_dir = Path(args.output_dir)

    try:
        adata = sc.read_h5ad(input_file)
    except Exception as e:
        print(f"Error reading file {input_file}: {e}")
        return

    print(f" Processing all samples together: {adata.n_obs} cells.")

    # Compute correlation
    try:
        compute_gene_correlation(
            adata=adata,
            correlation_type=args.correlation_type,
            replace_nans=args.replace_nans,
            min_cells_threshold=args.min_cells_threshold,
            file_path=input_file,
            output_dir=output_dir,
            xCellxSubject=args.xCellxSubject
        )
    except Exception as e:
        print(f" Error processing {input_file.name}): {e}")

    # Cleanup
    del adata
    gc.collect()


if __name__ == "__main__":
    main()


# compute_gene_correlation( adata,
# correlation_type='pearson', 
# replace_nans=True, 
# min_cells_threshold=20, 
# file_path="", 
# output_dir="",
# xCellxSubject=False)
