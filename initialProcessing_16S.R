##----------------------------------------------
## This script does the initial processing on the 
## BetterAir 16S data according to the DADA2 
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
setwd('~/Dropbox (BioBE)/BioBE Team Folder/Projects/BetterAir')


set.seed(2)
options(scipen=7)  # curtail scientific notation
options(digits=5) # number of digits to print on output

# install DADA2 for paired-end read assembly (instead of QIIME/Keaton's pipeline)
source("https://bioconductor.org/biocLite.R")
biocLite("dada2")
# also install ShortRead for DADA2 workflow
biocLite("ShortRead")
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
path <- "/Users/patrickfinnhorve/Desktop/16s_Dataset"

# list the files, sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1_001", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001", full.names = TRUE))

# check that there are equal numbers of files for forward and reverse
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

# extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_16S"), `[`, 1)

## examine quality profiles (random sample)
plotQualityProfile(fnFs[100:102]) + ylim(0,40) + geom_vline(xintercept=225)
plotQualityProfile(fnRs[12:14]) + ylim(0,40) + geom_vline(xintercept=160)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


##---------------------------
## Filtering and trimming: set truncLen based on where quality appeared to drop off, but leave enough for overlap=20
## takes ~2 minutes on 105 samples
##---------------------------

# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(225, 160),
#                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#                      compress=TRUE, multithread=TRUE) 
# head(out)

# filter and trim just the forward reads 
out <- filterAndTrim(fnFs, filtFs, truncLen=225,
                     trimLeft=15,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)


##---------------------------
## Learn error rates
##---------------------------

# Learn forward error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


##---------------------------
## Dereplicate sequences
##---------------------------
# NOTE: the KIT control had 0 reads after filtering/trimming, so it was removed from the list
filtFs <- filtFs[filtFs!="/Users/gwynhwyfer/Documents/ESBL/projects/BetterAir/BetterAir_rawData/16S/filtered/blue-treatment_KIT_F_filt.fastq.gz"]   
filtRs <- filtRs[filtRs!="/Users/gwynhwyfer/Documents/ESBL/projects/BetterAir/BetterAir_rawData/16S/filtered/blue-treatment_KIT_R_filt.fastq.gz"]   

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
sample.names <- sample.names[sample.names!="blue-treatment_KIT"]
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
## Merge paired reads - or don't, since they don't overlap enough (we used 319F - 806R primers)
##---------------------------

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


##---------------------------
## Construct sequence table
##---------------------------

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##---------------------------
## Remove chimeras
##---------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, file="seqtab.nochim_16S_noOverlap.rds")

##---------------------------
## Track reads through the pipeline
##---------------------------

#getN <- function(x) sum(getUniques(x))
#out <- out[row.names(out)!="blue-treatment_KIT_16S_S102_L001_R1_001.fastq.gz",]
#track <- data.frame("input"=out, "filtered"=sapply(dadaFs, getN), "denoisedF"=sapply(dadaRs, getN), "nonchim"=rowSums(seqtab.nochim))
#track$percentKept <- track[,5]/track[,1]
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "nonchim", "percentKept") #, "merged", "nonchim")
#rownames(track) <- sample.names
#head(track)


##---------------------------
## Assign taxonomy
##---------------------------
#seqtab.nochim <- readRDS(file="seqtab.nochim_noOverlap.rds")
taxa <- assignTaxonomy(seqtab.nochim, "/Users/gwynhwyfer/Documents/ESBL/projects/BetterAir/BetterAir_rawData/16S/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/gwynhwyfer/Documents/ESBL/projects/BetterAir/BetterAir_rawData/16S/tax/silva_species_assignment_v132.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(taxa, file="16Staxa_noOverlap.rds")

##---------------------------
## End pre-processing
##---------------------------

