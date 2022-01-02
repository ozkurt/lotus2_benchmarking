#Define the input (INP) and output (OUTP) dirs
set.dir(input=$INP)
set.dir(output=$OUTP)

######### DEMULTIPLEXING ###################
make.contigs(ffastq=$INP/A88BF.1.fq.gz, rfastq=$INP/A88BF.2.fq.gz,findex=$INP/A88BF.mid.fq.gz, rindex=$INP/A88BF.mid.fq.gz, oligos=oligos_A888.txt,bdiffs=2,pdiffs=2,processors=12)
make.contigs(ffastq=$INP/A739F.1.fq.gz, rfastq=$INP/A739F.2.fq.gz,findex=$INP/A739F.mid.fq.gz, rindex=$INP/A739F.mid.fq.gz, oligos=oligos_A739.txt,bdiffs=2,pdiffs=2,processors=12)


######### DENOISING & TAXONOMY ###################

#cat $OUTP/A88BF.1.fq.trim.contigs.fasta $OUTP/A739F.1.fq.trim.contigs.fasta > $OUTP/GUT.trim.contigs.fasta
#cat $OUTP/A88BF.1.fq.contigs.groups $OUTP/A739F.1.fq.contigs.groups > $OUTP/GUT.trim.contigs.groups

summary.seqs(fasta=$OUTP/GUT.trim.contigs.fasta)
screen.seqs(fasta=$OUTP/GUT.trim.contigs.fasta, group=$OUTP/GUT.trim.contigs.groups, maxambig=0, maxlength=301,processors=12)
summary.seqs(fasta=$OUTP/GUT.trim.contigs.good.fasta)

#Truncate the reads at 200 bases:
chop.seqs(fasta=$OUTP/GUT.trim.contigs.good.fasta, numbases=200, processors=12,keep=F)
unique.seqs(fasta=$OUTP/GUT.trim.contigs.good.chop.fasta)
count.seqs(name=$OUTP/GUT.trim.contigs.good.chop.names, group=$OUTP/GUT.trim.contigs.good.groups)
summary.seqs(count=$OUTP/GUT.trim.contigs.good.chop.count_table)

#Alignment:
align.seqs(fasta=$OUTP/GUT.trim.contigs.good.chop.unique.fasta, reference=silva.138.1.bacteria.gold.align)
summary.seqs(fasta=$OUTP/GUT.trim.contigs.good.chop.unique.align, count=$OUTP/GUT.trim.contigs.good.chop.count_table)



#Start and end positions are decided based on the summary file:
screen.seqs(fasta=$OUTP/GUT.trim.contigs.good.chop.unique.align, count=$OUTP/GUT.trim.contigs.good.chop.count_table, summary=$OUTP/GUT.trim.contigs.good.chop.unique.summary, start=15973, end=23444, processors=12)
summary.seqs(fasta=current, count=current)


#Filtering:
filter.seqs(fasta=$OUTP/GUT.trim.contigs.good.chop.unique.good.align, vertical=T, trump=.,processors=12)
unique.seqs(fasta=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.fasta, count=$OUTP/GUT.trim.contigs.good.chop.good.count_table)


pre.cluster(fasta=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.fasta, count=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.count_table, diffs=2)

#Chimera detection:
chimera.vsearch(fasta=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.fasta, count=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.count_table, dereplicate=t)
remove.seqs(fasta=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.fasta, accnos=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
summary.seqs(fasta=current, count=current)


#Classifying the sequences:
classify.seqs(fasta=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.pick.fasta, count=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax, cutoff=80)

#Removing the unknown taxa
remove.lineage(fasta=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.pick.fasta, count=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.pick.nr_v138_1.wang.taxonomy, taxon=unknown;)

summary.tax(taxonomy=current, count=current)

#OTU clustering:
dist.seqs(fasta=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03, processors=12)
cluster(column=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.pick.pick.dist, count=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
make.shared(list=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

classify.otu(list=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=$OUTP/GUT.trim.contigs.good.chop.unique.good.filter.unique.precluster.pick.nr_v138_1.wang.taxonomy, cutoff=80, label=0.03)

