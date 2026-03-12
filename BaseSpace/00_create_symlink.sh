#!/bin/bash

# Create symbolic links to FSTQ files

# BaseSpace requries files to have specific names, i.e. SampleName_S1_L001_R1_001.fastq.gz
# Create symbolic links with these filenames of our sequencing data

# Link arguments
LINK_ARG=files/symlink_args.txt

while read line
do
	ORIGINAL_PATH=`echo $line | cut -f 1 -d ' '`
	LINK_PATH=`echo $line | cut -f 2 -d ' '`

	ln -s $ORIGINAL_PATH $LINK_PATH

done < $LINK_ARG