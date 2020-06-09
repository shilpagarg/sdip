#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos='http://cran.us.r-project.org')
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse",repos='http://cran.us.r-project.org')
if (!requireNamespace("regioneR", quietly = TRUE))
  BiocManager::install("regioneR")
if (!requireNamespace("RSVSim", quietly = TRUE))
  BiocManager::install("RSVSim")

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(regioneR))
suppressPackageStartupMessages(library(RSVSim))

option_list = list(
  make_option(c("-d", "--dimensions"), action="store", type='character', help=".tsv file with chromosomes dimensions (as result from 'cut -f1,2 genome.fa.fai') [required]"),
  make_option(c("-n", "--numbers"), action="store", type='character', help="variant numbers for hom. DEL, het. DEL, hom. INS, het. INS, hom. SNP, het. SNP (e.g. -r '250:250:250:250:10000:10000')) [required]"),
  make_option(c("--hap1"), action="store", type='character', help="Output .tsv file with simulated SVs of haplotype 1 [required]"),
  make_option(c("--hap2"), action="store", type='character', help="Output .tsv file with simulated SVs of haplotype 2 [required]")
)


opt = parse_args(OptionParser(option_list=option_list))

if (is.null(opt$dimensions)) {
  stop('-d/--dimensions .tsv file is required')
} else {
  genome<-read.table(file.path(opt$dimensions), sep='\t', header = F)
}

variants = c('HOM DEL', 'HET DEL', 'HOM INS', 'HET INS', 'HOM SNP', 'HET SNP')

if (is.null(opt$numbers)) {
  stop('variant numbers in -n/--numbers are required')
}

numbers<-as.numeric(unlist(strsplit(opt$numbers, ':')))

if (length(variants) != length(numbers)) {
  stop('for each variant and genotype a ratio must be specified')
}

varwithproportions<-rep(variants, numbers)

del_lengths <- read.delim("/project/pacbiosv/data/simulated/simulated_data/scripts/simulation/del.lengths.txt", header=F)$V1*-1
del_lengths <- del_lengths[del_lengths >= 50 & del_lengths <= 100000]

hom_del_number <- length(which(varwithproportions == "HOM DEL"))
het_del_number <- length(which(varwithproportions == "HET DEL"))
hom_del_sizes = estimateSVSizes(n=hom_del_number, svSizes=del_lengths, minSize=50, maxSize=100000, hist=TRUE)
het_del_sizes = estimateSVSizes(n=het_del_number, svSizes=del_lengths, minSize=50, maxSize=100000, hist=TRUE)

ins_lengths <- read.delim("/project/pacbiosv/data/simulated/simulated_data/scripts/simulation/ins.lengths.txt", header=F)$V1
ins_lengths <- ins_lengths[ins_lengths >= 50 & ins_lengths <= 100000]

hom_ins_number <- length(which(varwithproportions == "HOM INS"))
het_ins_number <- length(which(varwithproportions == "HET INS"))
hom_ins_sizes = estimateSVSizes(n=hom_ins_number, svSizes=ins_lengths, minSize=50, maxSize=100000, hist=TRUE)
het_ins_sizes = estimateSVSizes(n=het_ins_number, svSizes=ins_lengths, minSize=50, maxSize=100000, hist=TRUE)

regions<-createRandomRegions(sum(numbers), length.mean=1, length.sd=0, genome=genome, non.overlapping=TRUE)

df_hap1 <- data.frame(chromosome=seqnames(regions),
                      start=start(regions),
                      end=end(regions)+1,
                      type=varwithproportions,
                      info=rep('INFO',sum(numbers)),
                      breakseqlen=rep(0, sum(numbers)),
                      genotype=rep('GT',sum(numbers)),
                      stringsAsFactors = FALSE)
df_hap2 <- data.frame(chromosome=c(),start=c(),end=c(),type=c(),info=c(),breakseqlen=c(),stringsAsFactors = FALSE)

num_hom_deletion <- 1
num_het_deletion <- 1
num_hom_insertion <- 1
num_het_insertion <- 1

for(i in (1:nrow(df_hap1))) {
  if (df_hap1$type[i] == 'HOM DEL') {
    df_hap1$end[i] <- df_hap1$start[i] + hom_del_sizes[num_hom_deletion]
    num_hom_deletion <- num_hom_deletion + 1
    df_hap1$info[i]<-'None'
    df_hap1$breakseqlen[i]<-sample(0:10,1)
    df_hap1$genotype[i]<-"1/1"
    df_hap1$type[i]<-"deletion"
    df_hap2 <- rbind(df_hap2, df_hap1[i,])
  } else if (df_hap1$type[i] == 'HET DEL') {
    df_hap1$end[i] <- df_hap1$start[i] + het_del_sizes[num_het_deletion]
    num_het_deletion <- num_het_deletion + 1
    df_hap1$info[i]<-'None'
    df_hap1$breakseqlen[i]<-sample(0:10,1)
    df_hap1$genotype[i]<-"1/0"
    df_hap1$type[i]<-"deletion"
  } else if (df_hap1$type[i] == 'HOM INS') {
    alphabet<-c('A', 'T', 'C', 'G')
    motif<-paste(sample(alphabet, hom_ins_sizes[num_hom_insertion], replace = T), collapse='')
    num_hom_insertion <- num_hom_insertion + 1
    df_hap1$info[i]<-motif
    df_hap1$breakseqlen[i]<-sample(0:10,1)
    df_hap1$genotype[i]<-"1/1"
    df_hap1$type[i]<-"insertion"
    df_hap2 <- rbind(df_hap2, df_hap1[i,])
  } else if (df_hap1$type[i] == 'HET INS') {
    alphabet<-c('A', 'T', 'C', 'G')
    motif<-paste(sample(alphabet, het_ins_sizes[num_het_insertion], replace = T), collapse='')
    num_het_insertion <- num_het_insertion + 1
    df_hap1$info[i]<-motif
    df_hap1$breakseqlen[i]<-sample(0:10,1)
    df_hap1$genotype[i]<-"1/0"
    df_hap1$type[i]<-"insertion"
  } else if (df_hap1$type[i] == 'HOM SNP') {
    alphabet<-c('A', 'T', 'C', 'G')
    motif<-paste(sample(alphabet, 1, replace = T), collapse='')
    df_hap1$info[i]<-motif
    df_hap1$genotype[i]<-"1/1"
    df_hap1$type[i]<-"SNP"
    df_hap2 <- rbind(df_hap2, df_hap1[i,])
  } else if (df_hap1$type[i] == 'HET SNP') {
     alphabet<-c('A', 'T', 'C', 'G')
     motif<-paste(sample(alphabet, 1, replace = T), collapse='')
     df_hap1$info[i]<-motif
     df_hap1$genotype[i]<-"1/0"
     df_hap1$type[i]<-"SNP"
   }
} 

options(scipen=999)
write.table(df_hap1,file = file.path(opt$hap1),quote = FALSE, col.names = FALSE, row.names = FALSE, sep='\t')
write.table(df_hap2,file = file.path(opt$hap2),quote = FALSE, col.names = FALSE, row.names = FALSE, sep='\t')


