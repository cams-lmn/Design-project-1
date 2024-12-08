import os
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.io import mmread

# Define the directory containing sample folders
samples_dir = 'processed_GEO_data'
output_stats_file = 'sample_stats.txt'
output_plots_dir = 'sample_plots'

# Ensure output directories exist
os.makedirs(output_plots_dir, exist_ok=True)

# Define QC thresholds for each sample
qc_thresholds = {
    'C1': {'nFeature_RNA': 2500, 'nCount_RNA': 5.0e3},
    'C2': {'nFeature_RNA': 2500, 'nCount_RNA': 5.0e3},
    'C3': {'nFeature_RNA': 1500, 'nCount_RNA': 5.0e3},
    'Early1': {'nFeature_RNA': 6000, 'nCount_RNA': 5.0e5},
    'Early2': {'nFeature_RNA': 3000, 'nCount_RNA': 5.0e4},
    'Early3': {'nFeature_RNA': 2500, 'nCount_RNA': 5.0e4}
}

# General QC threshold for all samples
general_threshold_nFeature_RNA = 200

# Function to load matrix data using scipy and pandas
def load_matrix_data(sample_dir, matrix_filename, features_filename, barcodes_filename):
    # Load the matrix (in Matrix Market format)
    matrix = mmread(os.path.join(sample_dir, matrix_filename)).tocsc()
    
    # Load features and barcodes (TSV files)
    features = pd.read_csv(os.path.join(sample_dir, features_filename), sep='\t', header=None)
    barcodes = pd.read_csv(os.path.join(sample_dir, barcodes_filename), sep='\t', header=None)
    
    return matrix, features, barcodes

# Initialize a dictionary to hold statistics
stats = {}

# Loop through each sample directory
for sample in os.listdir(samples_dir):
    sample_dir = os.path.join(samples_dir, sample)
    if os.path.isdir(sample_dir):
        print(f"Processing {sample}...")

        # Load RNA data
        rna_matrix, rna_features, barcodes = load_matrix_data(
            sample_dir, 'rna_matrix.mtx.gz', 'rna_features.tsv.gz', 'barcodes.tsv.gz')

        # Load ATAC data
        atac_matrix, atac_features, barcodes = load_matrix_data(
            sample_dir, 'atac_matrix.mtx.gz', 'atac_features.tsv.gz', 'barcodes.tsv.gz')

        # Compute statistics
        n_cells = barcodes.shape[0]
        n_rna_features = rna_features.shape[0]
        n_atac_features = atac_features.shape[0]

        # Sum the counts for RNA and ATAC
        rna_counts = np.array(rna_matrix.sum(axis=0)).flatten()
        atac_counts = np.array(atac_matrix.sum(axis=0)).flatten()

        # Calculate average counts per cell
        avg_rna_counts_per_cell = np.mean(rna_counts)
        avg_atac_counts_per_cell = np.mean(atac_counts)

        # Calculate the number of features (genes) detected per cell and average number of genes per cell
        features_per_cell_rna = np.array((rna_matrix > 0).sum(axis=0)).flatten()
        mean_features_per_cell_rna = np.mean(features_per_cell_rna)

        # Calculate the number of features (peaks) detected per cell and average number of peaks per cell
        features_per_cell_atac = np.array((atac_matrix > 0).sum(axis=0)).flatten()
        mean_features_per_cell_atac = np.mean(features_per_cell_atac)

        stats[sample] = {
            'Number of cells': n_cells,
            'Number of RNA features': n_rna_features,
            'Number of ATAC features': n_atac_features,
            'Average RNA UMI counts per cell': avg_rna_counts_per_cell,
            'Average ATAC counts per cell': avg_atac_counts_per_cell,
            'Average RNA features per cell': mean_features_per_cell_rna,
            'Average ATAC features per cell': mean_features_per_cell_atac
        }

        # Plot settings
        sns.set(style='whitegrid')
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        fig.suptitle(f'Distributions for {sample}')

        # RNA Counts per Cell
        sns.histplot(rna_counts, bins=50, kde=False, color="skyblue", ax=axes[0, 0], log_scale=(True, False))
        axes[0, 0].set_xlabel("RNA Counts per Cell (log scale)")
        axes[0, 0].set_ylabel("Frequency")
        axes[0, 0].set_title("Distribution of RNA Counts per Cell")
        axes[0, 0].axvline(qc_thresholds.get(sample, {}).get('nCount_RNA', np.nan), color='black', linestyle='-', linewidth=2, label=f"QC Threshold: {qc_thresholds.get(sample, {}).get('nCount_RNA', np.nan)}")
        axes[0, 0].axvline(avg_rna_counts_per_cell, color='red', linestyle='-', linewidth=2, label=f"Mean: {avg_rna_counts_per_cell:.2f}")
        axes[0, 0].legend()

        # ATAC Counts per Cell
        sns.histplot(atac_counts, bins=50, kde=False, color="orange", ax=axes[0, 1], log_scale=(True, False))
        axes[0, 1].set_xlabel("Accessibility Counts per Cell (log scale)")
        axes[0, 1].set_ylabel("Frequency")
        axes[0, 1].set_title("Distribution of ATAC Accessibility Counts per Cell")
        axes[0, 1].axvline(avg_atac_counts_per_cell, color='red', linestyle='-', linewidth=2, label=f"Mean: {avg_atac_counts_per_cell:.2f}")
        axes[0, 1].legend()

        # RNA Features per Cell
        sns.histplot(features_per_cell_rna, bins=50, kde=False, color="lightgreen", ax=axes[1, 0], log_scale=(True, False))
        axes[1, 0].set_xlabel("Number of Genes per Cell (log scale)")
        axes[1, 0].set_ylabel("Frequency")
        axes[1, 0].set_title("Distribution of Genes Detected per Cell (RNA)")
        axes[1, 0].axvline(general_threshold_nFeature_RNA, color='black', linestyle='-', linewidth=2, label='QC Threshold: 200')
        axes[1, 0].axvline(qc_thresholds.get(sample, {}).get('nFeature_RNA', np.nan), color='black', linestyle='-', linewidth=2, label=f"QC Threshold: {qc_thresholds.get(sample, {}).get('nFeature_RNA', np.nan)}")
        axes[1, 0].axvline(mean_features_per_cell_rna, color='red', linestyle='-', linewidth=2, label=f"Mean: {mean_features_per_cell_rna:.2f}")
        axes[1, 0].legend()

        # ATAC Features per Cell
        sns.histplot(features_per_cell_atac, bins=50, kde=False, color="purple", ax=axes[1, 1], log_scale=(True, False))
        axes[1, 1].set_xlabel("Number of Peaks per Cell (log scale)")
        axes[1, 1].set_ylabel("Frequency")
        axes[1, 1].set_title("Distribution of Peaks Detected per Cell (ATAC)")
        axes[1, 1].axvline(mean_features_per_cell_atac, color='red', linestyle='-', linewidth=2, label=f"Mean: {mean_features_per_cell_atac:.2f}")
        axes[1, 1].legend()
        
        # Save the plot with high quality
        plot_path = os.path.join(output_plots_dir, f"{sample}_distributions.png")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300)
        plt.close(fig)

# Save all statistics to output file
with open(output_stats_file, 'w') as f:
    for sample_name, sample_stats in stats.items():
        f.write(f"{sample_name}:\n")
        for stat_name, stat_value in sample_stats.items():
            f.write(f"  {stat_name}: {stat_value}\n")
        f.write("\n")

print("Analysis complete. Statistics and plots saved.")
