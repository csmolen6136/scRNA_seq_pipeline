#!/bin/bash

# Download DRAGEN single cell RNA output files from BaseSpace
# This script downloads the feature matrix, barcodes, and features (needed to interpret the matrix)

# NOTE: Ahead of uploading data, you will need to authenticate BaseSpace using $BASESPACE auth

# Path to BaseSpace CLI
BASESPACE=/data7/software/BaseSpace_Sequence_Hub/bs

# Get Appilcation ID
# You can also find Application IDs using the command: $BASESPACE appsession list
APP_ID=904474722

# Choose output directory
OUTPUT_DIR=/data7/organoid/scrna/data/DRAGEN_outputs/2026_06_03

$BASESPACE appsession download -i $APP_ID --exclude '*' --include '*.matrix.mtx.gz' --include '*.features.tsv.gz' --include '*.barcodes.tsv.gz' -o $OUTPUT_DIR --log=logs/2_download_outputs.log

