import argparse

import scanpy as sc
import scib

import matplotlib.pyplot as plt

# Use scanpy_env

# Perform clustering of scRNA-seq data

def clustering():
	# Parse arguments
	parser = argparse.ArgumentParser(description="scRNA-seq preprocessing")
	parser.add_argument("-o", "--output_dir", type=str, help="Output directory")
	args = parser.parse_args()

	OUTPUT_DIR=args.output_dir

	PROCESSED_DATA=f"{OUTPUT_DIR}/Data/Preprocessed_data.h5ad"
	OUTPUT_DATA=f"{OUTPUT_DIR}/Data/Uncorrected.h5ad"

	FIG_DIR=f"{OUTPUT_DIR}/Figures/Leiden_uncorrected"

	# Load pre-processed data
	adata=sc.read(PROCESSED_DATA)

	# Perform leiden clustering
	# Cluster over multiple resolutions to assess quality of clustering
	res_states=[0.2, 0.6, 0.8, 1, 1.2, 1.4, 1.6]
	for rs in res_states:
		sc.tl.leiden(adata, resolution=rs, key_added=f'leiden_{rs}')

	# Calculate silhouette scores
	sil_vals=[]
	for rs in res_states:
		sil=scib.me.silhouette(adata, label_key=f'leiden_{rs}', embed='X_pca')
		sil_vals.append(sil)
		sc.pl.umap(adata, color=f'leiden_{rs}',
					title=f'Leiden resolution {rs}: {sil}',
					frameon=False, legend_loc='on data',
					show=False)
		plt.savefig(f'{FIG_DIR}/Resolution_{rs}.png')
		plt.close()

	# Choose the resolution with the highest silhouette score to use for batch correction
	max_sil=max(sil_vals)
	max_idx=sil_vals.index(max_sil)
	reso=res_states[max_idx]
	print(f'Chose resolution {reso} with silhouette score {max_sil}')

	adata.obs["leiden"] = adata.obs[f"leiden_{reso}"]
	adata.obsm["Uncorrected"] = adata.obsm["X_pca"]

	# Save un-corrected data
	adata.write(OUTPUT_DATA, compression="gzip")

if __name__=="__main__":
	clustering()
