
library(dada2)
ncores=12
seqtab.nochim=list()
# Define a path for the dir to the 1st set of samples (path1)
# Define a path for the dir to the 2nd set of samples (path2)

path_vec=c(path1,path2)

for (i in 1:2){
path=path_vec[i]
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs = sort(list.files(path, pattern="1.fq.gz", full.names = TRUE))
fnRs = sort(list.files(path, pattern="2.fq.gz", full.names = TRUE))
# Extract sample names:
sample.names=sapply(strsplit(basename(fnFs), split='[.]'),`[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs = file.path(path, "filtered", paste0(sample.names, ".1.filt.fq.gz"))
filtRs = file.path(path, "filtered", paste0(sample.names, ".2.filt.fq.gz"))
names(filtFs) = sample.names
names(filtRs) = sample.names

# Since primer sequences cannot be removed in DADA2, they are trimmed from the reads of the soil-16S and the mock dataset using trimLeft = c(19, 34) and trimLeft=c(19,20) parameter, respectively. Otherwise:
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,truncQ=2, truncLen=200, rm.phix=TRUE,compress=TRUE, verbose=TRUE, multithread=ncores)

			  
#Learn the error rates:
errF <- learnErrors(filtFs, multithread=ncores)
errR <- learnErrors(filtRs, multithread=ncores)


#Sample Inference:
dadaFs <- dada(filtFs, err=errF, multithread=ncores,verbose=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=ncores,verbose=TRUE)

#Merge paired reads:
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#Construct sequence table:
seqtab <- makeSequenceTable(mergers)
seqtab.nochim[[i]] <- removeBimeraDenovo(seqtab, method="consensus", multithread= ncores, verbose=TRUE)

}


seqtab.nochim.merged=mergeSequenceTables(seqtab.nochim[[1]],seqtab.nochim[[2]])
#Chimera removal:
saveRDS(seqtab.nochim.merged, "seqtab.nochim_200.rds")
 
