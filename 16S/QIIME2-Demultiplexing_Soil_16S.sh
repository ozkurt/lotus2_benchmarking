#!/bin/bash

#Define the path for the already demultiplexed reads ($INP):

echo "Importing the already demultiplexed dataset:"
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $INP --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-soil.qza


echo "Removing primers:"

qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux-soil.qza \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r GGCCGYCAATTYMTTTRAGTTT \
--p-cores 12 \
--o-trimmed-sequences $INP/trimmed-seqs.qza \
--verbose

##### Only for Deblur: #######
echo "Joining the paired-end reads:"
qiime vsearch join-pairs --i-demultiplexed-seqs $INP/trimmed-seqs.qza --o-joined-sequences demux-soil-joined.qza
