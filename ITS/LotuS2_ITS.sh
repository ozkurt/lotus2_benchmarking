#!/bin/bash
#Define the directories for the map ($MAP), input ($INP) and output ($OUTP)
methods="cd-hit uparse vsearch"
sdm_params="sdm_miSeq_ITS.txt"



for method in $methods; do
	lotus2 -i $INP -o ${OUTP}_${method} -m $MAP -s $sdm_params -taxAligner lambda -refDB UNITE -amplicon_type ITS -tax_group fungi -ForwardPrimer ggGTGARTCATCRARTYTTTG -reversePrimer CCTSCSCTTANTDATATGC -CL $method -thr 12
done



