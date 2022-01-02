#!/bin/bash
#Define the directories for the map ($MAP), input ($INP) and output ($OUTP)

#For the soil-16S dataset, primers can be defined as a column in the mapping file or in the command: -forwardPrimer GTGYCAGCMGCCGCGGTAA -reversePrimer GGCCGYCAATTYMTTTRAGTTT
#-lulu was deactivated (-lulu 0) for the mock dataset.

methods="dada2 unoise uparse vsearch"
sdm_params="sdm_miSeq.txt"



for method in $methods; do
        for param_path in $sdm_params; do
            param_name="${param_path%.*}"
            lotus2 -i $INP -o ${OUTP}_${param_name}_${method} -m $MAP -s $param_path -taxAligner lambda -CL $method -thr 12
        done
done        
