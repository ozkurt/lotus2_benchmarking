#!/bin/bash

#Define an output dir (OUTP)

#"Demultiplexing" of the Soil-16S dataset slightly differs from the Gut-16S dataset. Please see QIIME2-Deblur_Soil_16S.sh for further information.
#"Denoising" and "Taxonomy" of the Gut- and Soil-16S datasets are performed using the same QIIME2 commands.


######### DEMULTIPLEXING ###################
# Import the dataset:
qiime tools import \
   --type EMPPairedEndSequences \
   --input-path A739_dataset \
   --output-path $OUTP/A739_dataset.qza


echo "Demultiplexing the A739 set:"


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


echo "Demultiplexing the A888 set"


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


echo "Joining the paired-end sequence reads:"

qiime vsearch join-pairs \
  --i-demultiplexed-seqs demux-A888.qza \
  --o-joined-sequences demux-A888-joined.qza


qiime vsearch join-pairs \
  --i-demultiplexed-seqs demux-A739.qza \
  --o-joined-sequences demux-A739-joined.qza


######### DENOISING ###################
echo "Starting to denoise with deblur:"

echo "Quality filtering:"

qiime quality-filter q-score \
 --i-demux demux-A888-joined.qza \
 --o-filtered-sequences $OUTP/demux-filtered-A888.qza \
 --o-filter-stats $OUTP/demux-filter-stats-A888.qza

echo "Denoising-deblur:"

qiime deblur denoise-16S --i-demultiplexed-seqs $OUTP/demux-filtered-A888.qza  --p-trim-length 200 --p-jobs-to-start 12 --o-representative-sequences $OUTP/rep-seqs-deblur-A888.qza --o-table $OUTP/table-deblur-A888.qza --p-sample-stats --o-stats $OUTP/deblur-stats-A888.qza


echo "Producing the visualization files:"


# Generate the sequence and OTU table:
qiime feature-table summarize \
  --i-table $OUTP/table-deblur-A888.qza \
  --o-visualization $OUTP/table-deblur-A888.qzv \
  --m-sample-metadata-file sampledata_A888.tsv

qiime feature-table tabulate-seqs \
  --i-data $OUTP/rep-seqs-deblur-A888.qza \
  --o-visualization $OUTP/rep-seqs-deblur-A888.qzv



qiime metadata tabulate \
  --m-input-file $OUTP/demux-filter-stats-A888.qza \
  --o-visualization $OUTP/demux-filter-stats-A888.qzv
qiime deblur visualize-stats \
  --i-deblur-stats $OUTP/deblur-stats-A888.qza \
  --o-visualization $OUTP/deblur-stats-A888.qzv


echo "Starting to denoise with deblur:"

echo "Quality filtering:"

qiime quality-filter q-score \
 --i-demux demux-A739-joined.qza \
 --o-filtered-sequences $OUTP/demux-filtered-A739.qza \
 --o-filter-stats $OUTP/demux-filter-stats-A739.qza

echo "Denoising-deblur:"

qiime deblur denoise-16S --i-demultiplexed-seqs $OUTP/demux-filtered-A739.qza  --p-trim-length 200 --p-jobs-to-start 12 --o-representative-sequences $OUTP/rep-seqs-deblur-A739.qza --o-table $OUTP/table-deblur-A739.qza --p-sample-stats --o-stats $OUTP/deblur-stats-A739.qza


echo "Producing the visualization files:"


# Generate the sequence and OTU table:
qiime feature-table summarize \
  --i-table $OUTP/table-deblur-A739.qza \
  --o-visualization $OUTP/table-deblur-A739.qzv \
  --m-sample-metadata-file sampledata_A739.tsv

qiime feature-table tabulate-seqs \
  --i-data $OUTP/rep-seqs-deblur-A739.qza \
  --o-visualization $OUTP/rep-seqs-deblur-A739.qzv



qiime metadata tabulate \
  --m-input-file $OUTP/demux-filter-stats-A739.qza \
  --o-visualization $OUTP/demux-filter-stats-A739.qzv
qiime deblur visualize-stats \
  --i-deblur-stats $OUTP/deblur-stats-A739.qza \
  --o-visualization $OUTP/deblur-stats-A739.qzv


echo "Mergeing tables.."

qiime feature-table merge \
        --i-tables $OUTP/table-deblur-A888.qza $OUTP/table-deblur-A739.qza \
        --o-merged-table $OUTP/table_deblur_merged.qza

echo "Mergeing seqs.."

qiime feature-table merge-seqs \
        --i-data $OUTP/rep-seqs-deblur-A888.qza $OUTP/rep-seqs-deblur-A739.qza \
        --o-merged-data $OUTP/rep-seqs-deblur_merged.qza


echo "Generating the outputs.."

# Generate the sequence and OTU table:
qiime feature-table summarize \
  --i-table $OUTP/table_deblur_merged.qza \
  --o-visualization $OUTP/table-merged-deblur.qzv \
  --m-sample-metadata-file sampledata_merged.tsv

qiime feature-table tabulate-seqs \
  --i-data $OUTP/rep-seqs-deblur_merged.qza \
  --o-visualization $OUTP/rep-seqs-merged-deblur.qzv


#Exporting the data
qiime tools export --input-path $OUTP/table_deblur_merged.qza --output-path $OUTP/exported-feature-tables

# Converting the biom file into txt file
biom convert -i $OUTP/exported-feature-tables/feature-table.biom -o $OUTP/table_merged_deblur.txt --to-tsv

######### TAXONOMY ###################
echo "Using already trained db"

echo "Classifying:"

qiime feature-classifier classify-sklearn \
  --i-classifier silva-nb-classifier.qza \
  --p-n-jobs 12\
  --i-reads $OUTP/rep-seqs-deblur_merged.qza  \
  --o-classification $OUTP/taxonomy_q2deblur.qza

echo "Importing the tax.qzv:"

qiime metadata tabulate \
  --m-input-file $OUTP/taxonomy_q2deblur.qza \
  --o-visualization $OUTP/taxonomy_q2deblur.qzv

