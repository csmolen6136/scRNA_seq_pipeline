import argparse

import pandas as pd
import numpy as np
import os

import scanpy as sc
import anndata as ad

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# Use scanpy_env

# Perform preprocessing and intial QC of scRNA-seq data

# Plot selected QC metrics
def plot_qc_metrics(anndata, pdf_name):
	pdf=PdfPages(pdf_name)

	# Genes and cells
	fig, axs=plt.subplots(ncols=2)
	sns.histplot(anndata.obs["n_genes_by_counts"], ax=axs[0], color='#8D6A9F')
	axs[0].set_title('Genes per cell')

	sns.histplot(anndata.var["n_cells_by_counts"], ax=axs[1], color='#B5F44A')
	axs[1].set_title('Cells per gene')
	pdf.savefig()

	# Mitochondrial count proportion
	fig, axs=plt.subplots(ncols=2)
	sns.histplot(anndata.obs["pct_counts_mito"], ax=axs[0], color='#F4B886')
	axs[0].set_title('Mitochondrial percentage')

	sns.histplot(anndata.obs["pct_counts_ribo"], ax=axs[1], color='#8CBCB9')
	axs[1].set_title('Ribosomal percentage')
	pdf.savefig()

	pdf.close()
	return

def preprocess():
	# Parse arguments
	parser = argparse.ArgumentParser(description="scRNA-seq preprocessing")
	parser.add_argument("-i", "--input_dir", type=str, help="Input directory")
	parser.add_argument("-o", "--output_dir", type=str, help="Output directory")
	parser.add_argument("--samp_data", type=str, help="Sample data table")
	parser.add_argument("-b", "--batch", type=str, help="Variable to use as batch")
	parser.add_argument("-v", "--var_cols", type=str, nargs='+',
					 		default=['Protocol', 'Extraction_batch', 'Sequencing_batch', 'Sample', 'Sex', 'Age'], help="Features to plot against UMAP")
	parser.add_argument("--min_genes", type=int, default=3000, help="Minimum genes in cell")
	parser.add_argument("--min_cells", type=int, default=10, help="Minimum cells expressing a gene")
	parser.add_argument("--mito", type=int, default=10, help="Maximum percentage of mitochondrial genes")
	parser.add_argument("--ribo", type=int, default=20, help="Maximum percentage of ribosomal genes")
	parser.add_argument("--var_genes", type=int, default=2000, help="Number of most highly variable genes to consider")
	args = parser.parse_args()

	INPUT_DIR=args.input_dir
	SAMPLE_INFO=args.samp_data
	OUTPUT_DIR=args.output_dir

	MIN_GENES=args.min_genes
	MIN_CELLS=args.min_cells
	MITO_PERC=args.mito
	RIBO_PERC=args.ribo
	N_HIGHLY_VARIABLE=args.var_genes

	BATCH_LABEL=args.batch
	rel_cols=args.var_cols

	# Define outputs
	QC_FIGURE=f'{OUTPUT_DIR}/Figures/QC_figure_prefilter.pdf'
	QC_FIGURE2=f'{OUTPUT_DIR}/Figures/QC_figure_postfilter.pdf'
	HV_FIG=f'{OUTPUT_DIR}/Figures/Highly_variable_genes.png'
	VAR_FIG=f'{OUTPUT_DIR}/Figures/PC_variance.pdf'
	UMAP_FIGURE=f'{OUTPUT_DIR}/Figures/UMAP_batches.png'

	COUNTS_TABLE=f'{OUTPUT_DIR}/Tables/Cell_counts.csv'

	PROCESSED_DATA=f'{OUTPUT_DIR}/Data/Preprocessed_data'

	# Load sample information
	samp_info=pd.read_csv(SAMPLE_INFO, sep='\t')

	# Find relevant directory for each sample
	for root, dirs, files in os.walk(INPUT_DIR):
		break
	samp_info['directory']=samp_info.SampleName.apply(lambda x: [i for i in dirs if f'{x}_ds' in i][0])

	# Load data for each sample
	adata={}
	for idx, row in samp_info.iterrows():
		dat=sc.read_10x_mtx(f'{INPUT_DIR}/{row.directory}',
							prefix=f'{row.SampleName}.scRNA.filtered.',
							make_unique=True) # Will warn about filename having multiple extensions - okay to ignore as long as it uses ".mtx" and ".gz"
		for rc in rel_cols:
			dat.obs[rc]=row[rc]

		adata[row.SampleName]=dat

	# Collapse data
	adata = ad.concat(adata, label="SampleName") # Will warn about non-unique names - we will fix in the next line
	adata.obs_names_make_unique()

	# Convert metadata to categorical
	for rc in rel_cols:
		adata.obs[rc]=adata.obs[rc].astype("category")

	# Save number of cells per sample and number of genes
	countdf=pd.DataFrame(adata.obs['SampleName'].value_counts())
	countdf.columns=['Starting_cells']

	# Annotate mitochondrial and ribosomal genes
	adata.var['mito']=adata.var_names.str.startswith("MT-")
	adata.var['ribo']=adata.var_names.str.startswith(("RPS", "RPL"))

	# Calculate QC metrics
	sc.pp.calculate_qc_metrics(adata, qc_vars=["mito", "ribo"], inplace=True, log1p=True)

	plot_qc_metrics(adata, QC_FIGURE)

	# Get percentage of cells per sample that fail thresholds
	mask=adata.obs['n_genes_by_counts']<MIN_GENES
	countdf['Low_expression']=countdf.index.map(adata.obs[mask].groupby('SampleName').size().to_dict())
	for i in range(2):
		vi=['mito', 'ribo'][i]
		thresh=[MITO_PERC, RIBO_PERC][i]
		col=['High_mito', 'High_ribo'][i]
		mask=adata.obs[f'pct_counts_{vi}']>thresh
		vi_count=adata.obs[mask].groupby('SampleName').size().to_dict()
		countdf[col]=countdf.index.map(vi_count)

	# Filter cells and genes and re-do QC plots
	# Filter on pre-defined thresholds
	sc.pp.filter_cells(adata, min_genes=MIN_GENES)
	sc.pp.filter_genes(adata, min_cells=MIN_CELLS)

	adata = adata[adata.obs['pct_counts_mito'] < MITO_PERC, :]
	adata = adata[adata.obs['pct_counts_ribo'] < RIBO_PERC, :]

	plot_qc_metrics(adata, QC_FIGURE2)

	# Save filtered cell counts
	countdf['Filtered_cells']=countdf.index.map(adata.obs['SampleName'].value_counts().to_dict())
	countdf.to_csv(COUNTS_TABLE)

	# Normalize
	adata.layers['counts']=adata.X.copy()
	sc.pp.normalize_total(adata) # Normalize counts per cell, such that all cells will have a count equal to the median of all cells
	sc.pp.log1p(adata) # Log + 1 normalize

	# Feature selection
	sc.pp.highly_variable_genes(adata, n_top_genes=N_HIGHLY_VARIABLE, batch_key=BATCH_LABEL)

	sc.pl.highly_variable_genes(adata)
	plt.savefig(HV_FIG)
	plt.close()

	# Dimensionality reduction
	sc.tl.pca(adata)

	# Plot PC variance
	sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
	plt.savefig(VAR_FIG)
	plt.close()

	# Compute nearest neighbors
	sc.pp.neighbors(adata, n_pcs=50)

	# UMAP and visualize
	sc.tl.umap(adata)

	sc.pl.umap(adata,
			color=['SampleName']+rel_cols,
			ncols=4)
	plt.savefig(UMAP_FIGURE)
	plt.close()

	# Save pre-processed data
	sc.write(PROCESSED_DATA, adata, compression='gzip')

if __name__=="__main__":
	preprocess()
	