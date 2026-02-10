#!/usr/bin/env python

import argparse
import json
from pathlib import Path
import scipy.sparse
import anndata as ad
import numpy as np
import gc

def slcGenes_and_cpm_adata(file_path, selected_genes, cpm_normalize=False, log_transform=False):
    """
    Filters AnnData to include only selected genes, adds missing genes with zero expression,
    and optionally applies CPM normalization and log transformation.
    """
    data_filtered = ad.read_h5ad(file_path)

    # Filter for selected genes present in the data
    current_genes = data_filtered.var_names.tolist()
    genes_to_keep = [gene for gene in selected_genes if gene in current_genes]
    data_filtered = data_filtered[:, genes_to_keep]

    # Identify and add missing genes
    missing_genes = [gene for gene in selected_genes if gene not in current_genes]

    if missing_genes:
        print(f"  Adding {len(missing_genes)} missing genes with zero expression")
        missing_expr = scipy.sparse.csr_matrix((data_filtered.shape[0], len(missing_genes)))
        data_filtered_new = ad.AnnData(
            X=scipy.sparse.hstack([missing_expr, data_filtered.X]).tocsc(),
            obs=data_filtered.obs.copy()
        )
        data_filtered_new.var_names = missing_genes + data_filtered.var_names.tolist()
        del data_filtered
        gc.collect()
    else:
        data_filtered_new = data_filtered

    data_filtered_new.X = scipy.sparse.csr_matrix(data_filtered_new.X)

    # CPM Normalization
    if cpm_normalize:
        print("  Performing CPM normalization")
        X = data_filtered_new.X
        if isinstance(X, scipy.sparse.spmatrix):
            total_counts = X.sum(axis=1).A1
        else:
            total_counts = X.sum(axis=1)
        total_counts[total_counts == 0] = 1
        data_filtered_new.X = (X.multiply(1e6 / total_counts[:, None])
                               if scipy.sparse.issparse(X) else X * 1e6 / total_counts[:, None])

    # Log Transformation
    if log_transform:
        print("  Performing log transformation")
        X = data_filtered_new.X
        if scipy.sparse.issparse(X):
            data_filtered_new.X = X.log1p()
        else:
            data_filtered_new.X = np.log1p(X)

    data_filtered_new.X = scipy.sparse.csr_matrix(data_filtered_new.X)
    return data_filtered_new

def main():
    parser = argparse.ArgumentParser(description='Select genes and apply CPM normalization')
    parser.add_argument('--input', required=True, help='Input h5ad file')
    parser.add_argument('--output', required=True, help='Output h5ad file')
    parser.add_argument('--gene_lists', required=True, help='Gene lists JSON file')
    parser.add_argument('--strategy', default='intersection', 
                        help='Gene selection strategy: intersection, union, or cell type name')
    parser.add_argument('--cpm_normalize', action='store_true', help='Apply CPM normalization')
    parser.add_argument('--no_log_transform', action='store_true', help='Skip log transformation')
    args = parser.parse_args()

    print(f"Processing file: {args.input}")

    # Load gene sets
    with open(args.gene_lists, "r") as fp:
        genes_set = json.load(fp)

    # Select gene list based on strategy
    if args.strategy == 'intersection':
        values = list(genes_set["less_stringent"].values())
        selected_genes = set(values[0])
        for value in values[1:]:
            selected_genes &= set(value)
        selected_genes = list(selected_genes)
        print(f"  Using intersection strategy: {len(selected_genes)} genes")
    elif args.strategy == 'union':
        selected_genes = list(set().union(*genes_set["less_stringent"].values()))
        print(f"  Using union strategy: {len(selected_genes)} genes")
    else:
        # Assume it's a cell type name
        if args.strategy in genes_set["less_stringent"]:
            selected_genes = genes_set["less_stringent"][args.strategy]
            print(f"  Using {args.strategy} gene list: {len(selected_genes)} genes")
        else:
            raise ValueError(f"Strategy '{args.strategy}' not found in gene lists")

    # Apply filtering and normalization
    filtered_data = slcGenes_and_cpm_adata(
        args.input,
        selected_genes,
        cpm_normalize=args.cpm_normalize,
        log_transform=not args.no_log_transform
    )

    # Save
    filtered_data.write_h5ad(args.output, compression="gzip")
    print(f"Saved processed file to: {args.output}")

if __name__ == "__main__":
    main()
