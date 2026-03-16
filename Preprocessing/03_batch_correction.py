import argparse

import scanpy as sc
import scanorama
import harmonypy as hm
import scvi

from torch.utils.data import DataLoader

from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

# Use scibm_env
# Needs to be run on laila

# Assess multiple batch correction methods

def batch_correction(adata, OUTPUT_DIR, BATCH_LABEL, SECOND=False):
	# Set output figures
	BENCHMARK_FIGURE=f"{OUTPUT_DIR}/Figures/scib_benchmark.pdf"
	BENCHMARK_TAB=f"{OUTPUT_DIR}/Tables/scib_benchmark.csv"
	suffix=''

	if SECOND:
		BENCHMARK_FIGURE=f"{OUTPUT_DIR}/Figures/scib_benchmark_2.pdf"
		BENCHMARK_TAB=f"{OUTPUT_DIR}/Tables/scib_benchmark_2.csv"
		adata.obsm["X_pca"]=adata.obsm["Batch_corrected"]
		suffix='_2'

	# Run multiple batch correction methods
	# Scanorama
	batch_cats = sorted(list(adata.obs[BATCH_LABEL].unique()))
	adata_list = [adata[adata.obs[BATCH_LABEL]==b].copy() for b in batch_cats]
	scanorama.integrate_scanpy(adata_list)

	adata.obsm[f"Scanorama{suffix}"] = np.zeros((adata.shape[0], adata_list[0].obsm["X_scanorama"].shape[1]))
	for i, b in enumerate(batch_cats):
		adata.obsm[f"Scanorama{suffix}"][adata.obs[BATCH_LABEL]  == b] = adata_list[i].obsm["X_scanorama"]

	# Harmony
	ho=hm.run_harmony(adata.obsm["X_pca"], adata.obs, BATCH_LABEL)
	adata.obsm[f"Harmony{suffix}"]=ho.Z_corr

	# scVI
	scvi.settings.dl_num_workers=6
	scvi.settings.num_threads=6
	scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=BATCH_LABEL)
	vae = scvi.model.SCVI(adata, gene_likelihood="nb", n_layers=2, n_latent=30)
	vae.train()
	adata.obsm[f"scVI{suffix}"] = vae.get_latent_representation()

	# scANVI
	lvae = scvi.model.SCANVI.from_scvi_model(
		vae,
		adata=adata,
		labels_key="leiden",
		unlabeled_category="Unknown",
	)
	lvae.train(max_epochs=20, n_samples_per_label=100)
	adata.obsm[f"scANVI{suffix}"] = lvae.get_latent_representation()

	# Benchmark
	bm = Benchmarker(
		adata,
		batch_key=BATCH_LABEL,
		label_key="leiden",
		bio_conservation_metrics=BioConservation(),
		batch_correction_metrics=BatchCorrection(),
		embedding_obsm_keys=["X_pca", f"Scanorama{suffix}", f"Harmony{suffix}", f"scVI{suffix}", f"scANVI{suffix}"],
		n_jobs=6)
	bm.benchmark()

	# Plot
	bm.plot_results_table()
	plt.savefig(BENCHMARK_FIGURE)
	plt.close()

	# Save to underlying dataframe
	df=bm.get_results(min_max_scale=False)
	df.to_csv(BENCHMARK_TAB)

	# Choose method with best overall score to continue downstream
	df=pd.read_csv(BENCHMARK_TAB)
	df=df[df.Embedding!='Metric Type']
	df.Total=df.Total.astype(float)
	best_meth=df[(df.Total==df.Total.max())].Embedding.to_list()[0]
	adata.obsm["Batch_corrected"]=adata.obsm[best_meth]

	# Save corrected data
	return adata

if __name__=="__main__":
	# Parse arguments
	parser = argparse.ArgumentParser(description="scRNA-seq preprocessing")
	parser.add_argument("-o", "--output_dir", type=str, help="Output directory")
	parser.add_argument("-b1", "--batch1", type=str, help="Variable to use as batch 1")
	parser.add_argument("-b2", "--batch2", type=str, default='', help="Variable to use as batch 2")
	args = parser.parse_args()

	OUTPUT_DIR=args.output_dir
	BATCH_LABEL1=args.batch1
	BATCH_LABEL2=args.batch2

	# Define inputs and outputs
	CLUSTERED_DATA=f"{OUTPUT_DIR}/Data/Uncorrected.h5ad"
	CORRECTED_DATA=f"{OUTPUT_DIR}/Data/Corrected.h5ad"

	# Load data
	adata=sc.read(CLUSTERED_DATA)

	# Correct on batch 1
	adata_corr=batch_correction(adata, OUTPUT_DIR, BATCH_LABEL1)

	# Correct on batch 2, if needed
	if BATCH_LABEL2!='':
		adata_corr=batch_correction(adata_corr, OUTPUT_DIR, BATCH_LABEL2, SECOND=True)

	# Save
	adata_corr.write(CORRECTED_DATA)
