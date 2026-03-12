import argparse

import scanpy as sc
import scib

import matplotlib.pyplot as plt

# Re-do clustering based on batch correction

def cluster_corrected():
	# Parse arguments
	parser = argparse.ArgumentParser(description="scRNA-seq preprocessing")
	parser.add_argument("-o", "--output_dir", type=str, help="Output directory")
	parser.add_argument("-v", "--var_cols", type=str, nargs='+',
					 		default=['Protocol', 'Extraction_batch', 'Sequencing_batch', 'Sample', 'Sex', 'Age'], help="Features to plot against UMAP")
	parser.add_argument("-bm", "--batch_method", type=str, default="", help="Alternative method to use for batch correction")
	args = parser.parse_args()

	OUTPUT_DIR=args.output_dir
	rel_cols=args.var_cols
	BM=args.batch_method

	CORRECTED_DATA=f"{OUTPUT_DIR}/Data/Corrected.h5ad"
	FIG_DIR=f"{OUTPUT_DIR}/Figures/Leiden_corrected"
	UMAP_FIGURE=f"{OUTPUT_DIR}/Figures/UMAP_batch_corrected.png"
	CLUSTERED_DATA=f"{OUTPUT_DIR}/Data/Clustered.h5ad"

	# Load data
	adata=sc.read_h5ad(CORRECTED_DATA)

	# Change batch correction method, if needed
	if BM!='':
		print('Batch correction method changed to', BM)
		adata.obsm["Batch_corrected"]=adata.obsm[BM]

	# Plot UMAP and batches
	# Compute nearest neighbors
	sc.pp.neighbors(adata, use_rep="Batch_corrected")

	# UMAP and visualize
	sc.tl.umap(adata)

	sc.pl.umap(adata,
			color=['SampleName']+rel_cols,
			ncols=4)
	plt.savefig(UMAP_FIGURE)
	plt.close()

	# After setting correction, re-do clustering
	res_states=[0.2, 0.6, 0.8, 1, 1.2, 1.4, 1.6]
	for rs in res_states:
		sc.tl.leiden(adata, resolution=rs, key_added=f'leiden_{rs}')

	# Calculate silhouette scores
	sil_vals=[]
	for rs in res_states:
		sil=scib.me.silhouette(adata, label_key=f'leiden_{rs}', embed='Batch_corrected')
		sil_vals.append(sil)
		sc.pl.umap(adata, color=f'leiden_{rs}',
					title=f'Leiden resolution {rs}: {sil}',
					frameon=False, legend_loc='on data',
					show=False)
		plt.savefig(f'{FIG_DIR}/Resolution_{rs}.png')
		plt.close()

	# Choose the resolution with the highest silhouette score to use going forward
	max_sil=max(sil_vals)
	max_idx=sil_vals.index(max_sil)
	reso=res_states[max_idx]
	print(f'Chose resolution {reso} with silhouette score {max_sil}')

	adata.obs["leiden"] = adata.obs[f"leiden_{reso}"]

	# Save clustered data
	adata.write(CLUSTERED_DATA)

if __name__=="__main__":
	cluster_corrected()

