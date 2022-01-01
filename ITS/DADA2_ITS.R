library(dada2)
ncores=12


# Define a path for the dir to the 1st set of samples (path1)
# Define a path for the dir to the 2nd set of samples (path2)

# Define a path for the UNITE db (unite.ref)


path_vec=c(path1, path2)
# First dataset: ####
for (i in 1:2) {
path=path_vec[i]
cutFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

filtFs <- file.path(path, "filtered", basename(cutFs))
filtRs <- file.path(path, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = ncores)
head(out)

#Learn the error rates:
errF <- learnErrors(filtFs, multithread = ncores)
errR <- learnErrors(filtRs, multithread = ncores)

#Dereplication:
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample Inference:
dadaFs <- dada(derepFs, err = errF, multithread = ncores)
dadaRs <- dada(derepRs, err = errR, multithread = ncores)

#Merge paired reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#Construct sequence table:
seqtab <- makeSequenceTable(mergers)
#Remove chimeras:
seqtab.nochim[[i]] <- removeBimeraDenovo(seqtab, method="consensus", multithread= ncores, verbose=TRUE)

}

seqtab_merged=mergeSequenceTables(seqtab.nochim[[1]], seqtab.nochim[[2]])
saveRDS(seqtab_merged,file="seqtab_merged_ITS.R")

taxa <- assignTaxonomy(seqtab_merged, unite.ref, multithread = ncores)
saveRDS(taxa,file="taxonomy_ITS_Rpkg_Specs.RDS")







