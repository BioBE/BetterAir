##----------------------------------------------
## This script does the initial processing on the 
## BetterAir ITS data according to the DADA2 
## statistical inference pipeline, including filtering  
## andmerging paired-end reads. 
## 
## Gwynne Mhuireach, 2018.10.05
## based on DADA2 tutorial for version 1.9.1
##----------------------------------------------


##---------------------------
## set up working directory, etc.
##---------------------------

setwd("/Users/gwynhwyfer/Documents/ESBL/projects/BetterAir")

set.seed(2)
options(scipen=7)  # curtail scientific notation
options(digits=5) # number of digits to print on output

# install DADA2 for paired-end read assembly (instead of QIIME/Keaton's pipeline)
source("https://bioconductor.org/biocLite.R")
biocLite("dada2")
# also install ShortRead for DADA2 workflow
biocLite("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
biocLite("DESeq2")
biocLite("phyloseq")

# The required package list:
reqpkg <- c("DESeq2", "ggplot2", "phyloseq", "lubridate", "zoo", "vegan", "ape",
            "xts", "VennDiagram", "stargazer", "dada2", "ShortRead")
# Load all required packages and show version
for (i in reqpkg) {
  print(i)
  print(packageVersion(i))
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}

# set ggplot2 theme elements
theme_set(theme_bw(base_size = 15))


##---------------------------
## file parsing for raw sequences (only need to do this if you haven't filtered sequences before)
##---------------------------

# CHANGE ME to the directory containing your demultiplexed forward-read fastq files (forward reads match *_R1_*)
path <- "/Users/gwynhwyfer/Documents/ESBL/projects/BetterAir/BetterAir_rawData/ITS"

# list the files, sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1_001", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001", full.names = TRUE))

# check that there are equal numbers of files for forward and reverse
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")


##---------------------------
## identify and remove primers
##---------------------------

FWD <- "CTTGGTCATTTAGAGGAAGTAA"  ## CHANGE ME to your forward primer sequence 
REV <- "GCTGCGTTCTTCATCGATGC"  ## CHANGE ME...

# verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# pre-filter sequences containing Ns
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# count number of times primers are found
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# install and use 'cutadapt' to remove primers
cutadapt <- "/Users/gwynhwyfer/.local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

  # count the presence of primers in the first cutadapt-ed sample
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


##---------------------------
## Filtering 
##---------------------------

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_ITS")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
# use this? sample.names <- sapply(strsplit(basename(fnFs), "_16S"), `[`, 1)
head(sample.names)

# visualize the quality profiles of the reads
plotQualityProfile(cutFs[1:4])
plotQualityProfile(cutRs[1:2])

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
#paste0(sample.names, "_F_filt.fastq.gz") use instead of basename(cutRs)

# Filter and trim
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 10), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)


##---------------------------
## Learn error rates
##---------------------------

# NOTE: the blue-treatment_CON1L and 'red-control_E2' samples had 0 reads after filtering/trimming, so they were removed from the list
filtFs <- filtFs[filtFs!="/Users/gwynhwyfer/Documents/ESBL/projects/BetterAir/BetterAir_rawData/ITS/cutadapt/filtered/blue-treatment_CON1L_ITS_S85_L001_R1_001.fastq.gz"]   
filtRs <- filtRs[filtRs!="/Users/gwynhwyfer/Documents/ESBL/projects/BetterAir/BetterAir_rawData/ITS/cutadapt/filtered/blue-treatment_CON1L_ITS_S85_L001_R2_001.fastq.gz"]   
filtFs <- filtFs[filtFs!="/Users/gwynhwyfer/Documents/ESBL/projects/BetterAir/BetterAir_rawData/ITS/cutadapt/filtered/red-control_E2_ITS_S48_L001_R1_001.fastq.gz"]   
filtRs <- filtRs[filtRs!="/Users/gwynhwyfer/Documents/ESBL/projects/BetterAir/BetterAir_rawData/ITS/cutadapt/filtered/red-control_E2_ITS_S48_L001_R2_001.fastq.gz"]   

# Learn forward error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


##---------------------------
## Dereplicate sequences
##---------------------------

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
sample.names <- sample.names[sample.names!="blue-treatment_CON1L"]
sample.names <- sample.names[sample.names!="red-control_E2"]
names(derepFs) <- sample.names
names(derepRs) <- sample.names


##---------------------------
## Sample inference
##---------------------------

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
length(derepFs)
dadaFs[[1]]
length(dadaFs)
length(dadaRs)

##---------------------------
## Don't merge paired reads, use forward only bc reverse were poor quality
##---------------------------

#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate = TRUE)
# Inspect the merger data.frame from the first sample
#head(mergers[[1]])


##---------------------------
## Construct sequence table
##---------------------------

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##---------------------------
## Remove chimeras
##---------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, file="seqtab.nochim_ITS_fReadsOnly.rds")

##---------------------------
## Track reads through the pipeline
##---------------------------

getN <- function(x) sum(getUniques(x))
length(row.names(out2))
out2 <- out[row.names(out)!="blue-treatment_CON1L_ITS_S85_L001_R1_001.fastq.gz",]
out2 <- out2[row.names(out2)!="red-control_E2_ITS_S48_L001_R1_001.fastq.gz",]

track <- data.frame("input"=out2, "filtered"=sapply(dadaFs, getN), "denoisedF"=sapply(dadaRs, getN), "nonchim"=rowSums(seqtab.nochim))
track$percentKept <- track[,5]/track[,1]
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "nonchim", "percentKept") #, "merged", "nonchim")
rownames(track) <- sample.names
head(track)


##---------------------------
## Assign taxonomy
##---------------------------
#seqtab.nochim <- readRDS(file="seqtab.nochim_ITS_fReads.rds")
unite.ref <- "/Users/gwynhwyfer/Documents/ESBL/projects/BetterAir/BetterAir_rawData/ITS/tax/sh_general_release_dynamic_s_01.12.2017.fasta"  
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(taxa, file="ITStaxa_fReadsOnly.rds")

##---------------------------
## End pre-processing
##---------------------------


