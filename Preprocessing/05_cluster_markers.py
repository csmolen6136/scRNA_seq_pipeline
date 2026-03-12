import argparse

import scanpy as sc
import matplotlib.pyplot as plt

# Identify marker genes for each cell cluster

def marker_genes():
	# Parse arguments
	parser = argparse.ArgumentParser(description="scRNA-seq preprocessing")
	parser.add_argument("-o", "--output_dir", type=str, help="Output directory")
	parser.add_argument("-r", "--resolution", type=str, default="", help="Alternative resolution to use for defined clusters")
	parser.add_argument("-b", "--batch", type=str, help="Variable to use as batch")
	parser.add_argument("-g", "--genes", type=str, nargs='+',
					 		default=['CDR2', 'POLR3E', 'EEF2K', 'SDR42E2', 'VWA3A', 'MOSMO', 'PDZD9', 'UQCRC2'], help="Specific genes to track across clusters")
	args = parser.parse_args()

	OUTPUT_DIR=args.output_dir
	RES=args.resolution
	BATCH=args.batch

	PLOT_GENES=args.genes

	CLUSTERED_DATA=f"{OUTPUT_DIR}/Data/Clustered.h5ad"

	GENE_LIST=f"{OUTPUT_DIR}/Tables/Cluster_marker_genes.csv"
	MARKER_DOTPLOT=f"{OUTPUT_DIR}/Figures/Cluster_markers.pdf"
	INPUT_GENE_PLOT=f"{OUTPUT_DIR}/Figures/Cluster_gene_list.png"

	# Load data
	adata=sc.read_h5ad(CLUSTERED_DATA)

	# Update resolution, if needed
	if RES!='':
		adata.obs["leiden"] = adata.obs[f"leiden_{RES}"]
	
	# Get cluster-specific differentially expressed genes
	sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
	markerdf=sc.get.rank_genes_groups_df(adata, group=None, log2fc_min=2)
	markerdf.sort_values(by=['group', 'pvals_adj', 'logfoldchanges'], ascending=[True, True, False], inplace=True)
	markerdf.to_csv(GENE_LIST, index=False)

	# Plot
	sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden", n_genes=5, min_logfoldchange=2, show=False)
	plt.savefig(MARKER_DOTPLOT, bbox_inches='tight')
	plt.close()

	# Plot expression of specific genes
	sc.pl.umap(adata,
				color=PLOT_GENES+['leiden', BATCH],
				legend_loc="on data",
    			frameon=False,
    			ncols=4, show=False)
	plt.tight_layout()
	plt.savefig(INPUT_GENE_PLOT)
	plt.close()

if __name__=="__main__":
	marker_genes()