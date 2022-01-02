#!/bin/bash

#Define an output dir (OUTP):


######### DEMULTIPLEXING ###################
# Import the dataset:
qiime tools import \
   --type EMPPairedEndSequences \
   --input-path A739_dataset \
   --output-path $OUTP/A739_dataset.qza


echo "Demultiplexing the A739 dataset:"


# Demultiplex the reads:
qiime demux emp-paired \
  --m-barcodes-file sampledata_A739.tsv \
  --m-barcodes-column barcode-sequence \
  --p-rev-comp-mapping-barcodes \
  --i-seqs $OUTP/A739_dataset.qza \
  --verbose \
  --o-per-sample-sequences $OUTP/demux-A739.qza \
  --o-error-correction-details $OUTP/demux-details-A739.qza




# Visualize the demultiplexed data:
qiime demux summarize \
  --i-data $OUTP/demux-A739.qza \
  --o-visualization $OUTP/demux-A739.qzv



# Import the dataset:
qiime tools import \
   --type EMPPairedEndSequences \
   --input-path A888_dataset \
   --output-path $OUTP/A888_dataset.qza


echo "Demultiplexing the A888 dataset:"


# Demultiplex the reads:
qiime demux emp-paired \
  --m-barcodes-file sampledata_A888.tsv \
  --m-barcodes-column barcode-sequence \
  --p-rev-comp-mapping-barcodes \
  --i-seqs $OUTP/A888_dataset.qza \
  --verbose \
  --o-per-sample-sequences $OUTP/demux-A888.qza \
  --o-error-correction-details demux-details-A888.qza


# Visualize the demultiplexed data:
qiime demux summarize \
  --i-data $OUTP/demux-A888.qza \
  --o-visualization $OUTP/demux-A888.qzv
  

######### DENOISING ###################
echo "Denoising with DADA2:"

#Denoising of the demultiplexed reads with dada2:

#Denoising the A739 dataset:
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-A739.qza \
  --p-trunc-len-f 200 \
  --p-trunc-len-r 200 \
  --p-n-threads 12 \
  --verbose \
  --o-table $OUTP/table-A739.qza \
  --o-representative-sequences $OUTP/rep-seqs-A739.qza \
  --o-denoising-stats $OUTP/denoising-stats-A739.qza


echo "Generating the outputs:"

# Generate the sequence and OTU table:
qiime feature-table summarize \
  --i-table $OUTP/table-A739.qza \
  --o-visualization $OUTP/table-A739.qzv \
  --m-sample-metadata-file sampledata_A739.tsv

qiime feature-table tabulate-seqs \
  --i-data $OUTP/rep-seqs-A739.qza \
  --o-visualization $OUTP/rep-seqs-A739.qzv


# Visualize the denoising stats:
qiime metadata tabulate \
  --m-input-file $OUTP/denoising-stats-A739.qza \
  --o-visualization $OUTP/denoising-stats-A739.qzv



# Denoising the A888 dataset:
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-A888.qza \
  --p-trunc-len-f 200 \
  --p-trunc-len-r 200 \
  --p-n-threads 12 \
  --verbose \
  --o-table $OUTP/table-A888.qza \
  --o-representative-sequences $OUTP/rep-seqs-A888.qza \
  --o-denoising-stats $OUTP/denoising-stats-A888.qza


echo "Generating the outputs:"

# Generate the sequence and ASV table:
qiime feature-table summarize \
  --i-table $OUTP/table-A888.qza \
  --o-visualization $OUTP/table-A888.qzv \
  --m-sample-metadata-file sampledata_A888.tsv

qiime feature-table tabulate-seqs \
  --i-data $OUTP/rep-seqs-A888.qza \
  --o-visualization $OUTP/rep-seqs-A888.qzv


# Visualize the denoising stats:
qiime metadata tabulate \
  --m-input-file $OUTP/denoising-stats-A888.qza \
  --o-visualization $OUTP/denoising-stats-A888.qzv

echo "Mergeing tables:"

qiime feature-table merge \
        --i-tables $OUTP/table-A739.qza $OUTP/table-A888.qza \
        --o-merged-table $OUTP/table_merged.qza

echo "Mergeing seqs:"

qiime feature-table merge-seqs \
        --i-data $OUTP/rep-seqs-A739.qza $OUTP/rep-seqs-A888.qza \
        --o-merged-data $OUTP/rep-seqs_merged.qza


echo "Generating the outputs:"

# Generate the sequence and OTU table:
qiime feature-table summarize \
  --i-table $OUTP/table_merged.qza \
  --o-visualization $OUTP/table-merged-dada2.qzv \
  --m-sample-metadata-file sampledata_merged.tsv

qiime feature-table tabulate-seqs \
  --i-data $OUTP/rep-seqs_merged.qza \
  --o-visualization $OUTP/rep-seqs-merged-dada2.qzv


#Exporting the data
qiime tools export --input-path $OUTP/table_merged.qza --output-path $OUTP/exported-feature-tables

# Converting the biom file into txt file
biom convert -i $OUTP/exported-feature-tables/feature-table.biom -o $OUTP/table_merged_d2.txt --to-tsv



######### TAXONOMY ###################
echo "Using an already trained db"

echo "Classifying:"

qiime feature-classifier classify-sklearn \
  --i-classifier silva-nb-classifier.qza \
  --p-n-jobs 12\
  --i-reads rep-seqs_merged.qza \
  --o-classification taxonomy_d2.qza

echo "Importing the tax.qzv:"

qiime metadata tabulate \
  --m-input-file taxonomy_d2.qza \
  --o-visualization taxonomy_d2.qzv

