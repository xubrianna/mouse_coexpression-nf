import argparse
import scanpy as sc
import anndata as ad
import gc
from pathlib import Path
import numpy as np
import scipy.sparse
import pandas as pd 
from scipy.stats import spearmanr
import bottleneck as bn
import psutil
import os 
import sys

def rank(data):
    """Rank normalize data with NaN handling"""
    orig_shape = data.shape
    data = bn.nanrankdata(data).reshape(-1) - 1  # Ranks, ignoring NaNs
    # Normalize between 0 and 1
    return (data / np.nansum(~np.isnan(data))).reshape(orig_shape).astype(np.float32)

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
        rank_normalized_matrix[np.isnan(rank_normalized_matrix)] = bn.nanmean(rank_normalized_matrix)

    print(f"Finished processing sample {sample_id}")
    del corr_matrix, subset_data
    gc.collect()
    return sample_id, {
        'rank_normalized_matrix': rank_normalized_matrix,
        'gene_names': adata.var_names.tolist()
    }

def write_edge_list_chunked(matrix, gene_names, output_file, value_col='rank_normalized_correlation',
                            upper_only=True, chunk_rows=50):
    """Write edge list to CSV in row-wise chunks to avoid large intermediate DataFrames.

    Args:
        matrix: 2-D numpy array (n_genes x n_genes)
        gene_names: sequence of gene name strings
        output_file: path-like or str
        value_col: column name for the correlation values
        upper_only: if True write upper triangle only; if False write all off-diagonal entries
        chunk_rows: number of source-gene rows to process per chunk
    """
    n = len(gene_names)
    gene_arr = np.asarray(gene_names, dtype='object')
    header_written = False

    for i_start in range(0, n, chunk_rows):
        i_end = min(i_start + chunk_rows, n)
        r_list, c_list, v_list = [], [], []
        
        for i in range(i_start, i_end):
            if upper_only:
                j_start = i + 1
                if j_start >= n:
                    continue
                j_arr = np.arange(j_start, n, dtype=np.int32)
            else:
                j_arr = np.concatenate([
                    np.arange(0, i, dtype=np.int32),
                    np.arange(i + 1, n, dtype=np.int32)
                ])
            if len(j_arr) == 0:
                continue
            r_list.append(np.full(len(j_arr), i, dtype=np.int32))
            c_list.append(j_arr)
            v_list.append(matrix[i, j_arr])

        if not r_list:
            continue

        r = np.concatenate(r_list)
        c = np.concatenate(c_list)
        v = np.concatenate(v_list)
        
        # Build DataFrame minimally and write immediately
        chunk_df = pd.DataFrame({
            'geneA': gene_arr[r],
            'geneB': gene_arr[c],
            value_col: v.astype(np.float32)
        })
        chunk_df.set_index(['geneA', 'geneB'], inplace=True)
        chunk_df.to_csv(output_file, mode='a', header=not header_written)
        header_written = True
        
        del r, c, v, chunk_df, r_list, c_list, v_list
        gc.collect()


def cor_mtx_to_edgeList(file_path, agg_adata, edge_list_dir2):
    """Write full edge list (all off-diagonal entries) with sanity check."""
    n = agg_adata.shape[0]
    # Sanity check: every gene should appear as source with n-1 targets (no NaN off-diagonals)
    matrix = agg_adata.X
    off_diag_counts = n - 1 - np.sum(np.isnan(matrix), axis=1)
    
    if not np.all(off_diag_counts == n - 1):
        print(f"Sanity check failed ({np.sum(off_diag_counts != n-1)} genes have NaNs). Skipping file.\n")
        return
    
    print("Sanity check passed")
    output_file = edge_list_dir2 / f"{file_path.stem}_edgeList.csv"
    
    # Write directly without intermediate DataFrame copies
    write_edge_list_chunked(matrix, agg_adata.var_names, output_file,
                            value_col='rank_norm_corr', upper_only=False)
    print(f"Saved edge list to: {output_file}\n")

def correlation_matrix_to_edgelist(matrix_np, gene_names, output_file):
    """Write upper-triangle edge list directly to CSV (memory-efficient)."""
    write_edge_list_chunked(matrix_np, gene_names, output_file,
                            value_col='rank_normalized_correlation', upper_only=True)

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

    #1. Save AnnData object (convert to float32 if not already)
    aggregated_matrix = cor["aggregated_rank_normalized_matrix"].astype(np.float32)

    gene_names = cor["gene_names"]
    agg_adata = ad.AnnData(X=aggregated_matrix)
    agg_adata.var_names = gene_names
    agg_adata.write_h5ad(agg_file, compression="gzip")
    print(f"Saved aggregated matrix and gene names to {agg_file}")

    #2. Save upper-triangle edgeList
    correlation_matrix_to_edgelist(aggregated_matrix, gene_names, edge_list_file)
    print("saved upper-triangle edge list")
    del aggregated_matrix  # Release aggregated matrix early after edge list writing
    gc.collect()

    #2.2 Save full edgeList (all off-diagonal)
    cor_mtx_to_edgeList(file_path, agg_adata, edge_list_dir2)
    del agg_adata  # Release AnnData object after edge list writing
    gc.collect()

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
        aggregated_matrix = rank(aggregated_matrix).astype(np.float32)
    else:
        print("No valid samples to aggregate!")
        return

    # Summary
    summary_df = pd.DataFrame({
        "total_genes": [len(gene_names)],
        "total_single_cells": [subset.shape[0]],
        "total_samples": [sample_count]
    })

    category_data = {
        'aggregated_rank_normalized_matrix': aggregated_matrix,
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

    # Clean up memory - note: aggregated_matrixEdge was never defined, removed
    del subset, aggregated_matrix, category_data, summary_df
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
        sys.exit(1)

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
        sys.exit(1)

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
