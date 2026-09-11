# setwd("~/Desktop/Region-analysis")
library("regioneR")
set.seed(12345)
A <- toGRanges("data/Bl-pirs-all-matches-hg18-6c-lifted-orig-hg38-2408-rem-random.bed") #(This should be replaced by your interest of RNA bed file)
B <- toGRanges("data/chr1-22.Gene-2000UP.bed") #(This should be replaced by your interest of genome regions)
#B <- filterChromosomes(B, chr.type="autosomal", organism="hg19")
pt <- overlapPermTest(A, B, ntimes=100, genome="hg38", count.once=TRUE, alternative="auto")
pt
summary(pt$numOverlaps)

#use the following:
#chr1-22.gene.bed
#chr1-22.intron.bed 
#chr1-22.codingexon.bed 
#chr1-22.exonplus.bed 
#chr1-22.CpG.bed 
#chr1-22.5utr.bed
#chr1-22.3utr.bed
#chr1-22.Gene-2000Down.bed
#chr1-22.Gene-2000UP.bed