#!/bin/bash

# Upload FASTQ files to BaseSpace

# NOTE: Ahead of uploading data, you will need to authenticate BaseSpace using $BASESPACE auth

# NOTE 2: DRAGEN requires filenames to fit specific formats. Filenames must be: SampleName_S1_L001_R1_001.fastq.gz
# SampleName cannot contain underscores ( _ )
# Script 00_create_symlink.sh creates symbolic links between the original FASTQ files (incorrect names) and the new filenames (corrected names)
# This is a way to provide the correct filename without renaming or copying the original file

# Path to BaseSpace CLI
BASESPACE=/data7/software/BaseSpace_Sequence_Hub/bs

# Path to FASTQ files (or symbolic links to FASTQ files)
FASTQ_DIR=sym_links/

# Get project ID
# Projects are set up in Illumina BaseSpace on cloud
# You can also find them using the command: $BASESPACE list project
PROJECT_ID=491500010

$BASESPACE dataset upload -p $PROJECT_ID --recursive $FASTQ_DIR --log=logs/1_upload_data.log 
