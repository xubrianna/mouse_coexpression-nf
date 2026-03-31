#!/usr/bin/env python

import argparse
import json
from pathlib import Path
from itertools import combinations
import anndata as ad

def analyze_gene_dict(gene_dict, frequency_threshold=0.2):
    """
    Processes a gene dictionary to generate:
    1. A list of genes found per cell type across all studies (stringent).
    2. A list of genes found per cell type in at least __ of the studies (less stringent).
    3. A dictionary of pairwise intersections of "less stringent" gene lists between cell types.
    
    Args:
        gene_dict (dict): The input gene dictionary {cell_type -> {study_name -> [genes]}}.
    
    Returns:
        dict: Results with keys:
            - 'stringent': Genes found across all studies per cell type.
            - 'less_stringent': Genes found in at least __ of studies per cell type.
            - 'pairwise_intersections': Intersections of "less stringent" lists between cell types.
    """
    results = {
        "stringent": {},
        "less_stringent": {},
        "pairwise_intersections": {}
    }

    # Process each cell type in the gene dictionary
    for cell_type, studies in gene_dict.items():
        study_genes = list(studies.values())
        num_studies = len(study_genes)

        # Compute stringent list: Genes found in all studies
        stringent_genes = set(study_genes[0])
        for genes in study_genes[1:]:
            stringent_genes.intersection_update(genes)
        results["stringent"][cell_type] = list(stringent_genes)

        # Compute less stringent list: Genes found in at least frequency_threshold of studies
        gene_counts = {}
        for genes in study_genes:
            for gene in genes:
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
        less_stringent_genes = [
            gene for gene, count in gene_counts.items()
            if count >= frequency_threshold * num_studies
        ]
        results["less_stringent"][cell_type] = less_stringent_genes

    # Compute pairwise intersections of less stringent gene lists between cell types
    cell_types = list(results["less_stringent"].keys())
    for cell_type1, cell_type2 in combinations(cell_types, 2):
        intersect_genes = set(results["less_stringent"][cell_type1]).intersection(
            results["less_stringent"][cell_type2]
        )
        results["pairwise_intersections"][(cell_type1, cell_type2)] = list(intersect_genes)

    return results

def extract_metadata_from_filename(filename):
    """
    Extract cell type and study name from a filename.
    Assumes the format: "CellType_StudyName_Year.h5ad".
    """
    name_parts = filename.split('_')
    cell_type = name_parts[0]
    study_name = name_parts[-1].split('.')[0]
    return cell_type, study_name

def process_h5ad_files(directory):
    """
    Process all .h5ad files in the specified directory, extract metadata,
    and populate a dictionary with cell type and study information.
    """
    directory = Path(directory)  # Ensure the directory is a Path object
    gene_dict = {} # Dictionary to store metadata information

    for file_path in directory.glob("*.h5ad"):
        print(f"\nProcessing file: {file_path.name}")
        # Load the AnnData object
        adata = ad.read_h5ad(file_path)
        print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes.")

        # Extract metadata (cell type, study name)      
        cell_type, study_name = extract_metadata_from_filename(file_path.name)
        # Populate the dictionary
        if cell_type not in gene_dict:
            gene_dict[cell_type] = {}
        if study_name not in gene_dict[cell_type]:
            gene_dict[cell_type][study_name] = []

        gene_dict[cell_type][study_name] = adata.var.index.tolist()
        print(f"Stored {len(adata.var.index)} genes for {cell_type} - {study_name}.")

    return gene_dict

# Utility function for JSON serialization
def convert_tuple_keys(obj):
    """Convert tuple keys to strings for JSON serialization."""
    if isinstance(obj, dict):
        return {str(k): convert_tuple_keys(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_tuple_keys(item) for item in obj]
    else:
        return obj

def main():
    parser = argparse.ArgumentParser(description='Generate gene lists from h5ad files')
    parser.add_argument('--input_dir', required=True, help='Input directory with h5ad files')
    parser.add_argument('--output', required=True, help='Output JSON file')
    parser.add_argument('--frequency_threshold', type=float, default=0.2, 
                        help='Frequency threshold for less stringent genes')
    args = parser.parse_args()

    print("\nExtracting genes...")
    gene_dict = process_h5ad_files(args.input_dir)

    print("\nAnalyzing gene intersections...")
    results = analyze_gene_dict(gene_dict, args.frequency_threshold)
    results_serializable = convert_tuple_keys(results)

    print("\nSaving results to JSON...")
    with open(args.output, 'w') as f:
        json.dump(results_serializable, f, indent=4)

    print("\nSummary of stringent genes per cell type:")
    for cell_type, genes in results["stringent"].items():
        print(f"{cell_type}: {len(genes)} genes")

    print("\nSummary of less stringent genes per cell type:")
    for cell_type, genes in results["less_stringent"].items():
        print(f"{cell_type}: {len(genes)} genes")

if __name__ == "__main__":
    main()
