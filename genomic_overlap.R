#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
library(GenomicRanges)
snps<- read.table(args[1],header=T)
gr <- makeGRangesFromDataFrame(snps, TRUE)
d <- disjoin(gr)
olaps <- findOverlaps(d, gr)
mcols(d) <- splitAsList(gr$sample[subjectHits(olaps)], queryHits(olaps))
out <- d[elementLengths(d$X) > 1]
df <- data.frame(chr=seqnames(out),segment=ranges(out))
write.table(df,file=args[2],quote=F,sep="\t",row.names=F)
