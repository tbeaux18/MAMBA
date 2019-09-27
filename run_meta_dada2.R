#!/usr/bin/env Rscript
library(tidyverse)
library(ggplot2)
library (magrittr)
library(stringr)
library(dada2) # Rcpp # only available for R 3.6 install via biocmanager 
library(phyloseq) # install via biocmanager 
library(Biostrings) # install via biocmanager 


# fastq.path <- "/data/MAMBA/meta_fastq" # CHANGE ME to the directory containing the fastq files after unzipping.
fastq.path <- "~/Documents/loyola/MAMBA"
list.files(fastq.path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(fastq.path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(fastq.path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# need to create a directory and paths to print these plots
# plots the quality scores
plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])

# filtFq.path <- "/data/MAMBA/filtered"
mamba.path <- "~/Documents/loyola/MAMBA"
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(mamba.path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(mamba.path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# figures out usuable reads
# need to determine aggregate where to truncate
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,230),
                     trimLeft=10, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
dout <- as.data.frame(out)
FaT.matrix <- file.path(mamba.path, "filtered_trimmed_reads.csv")
write.csv(dout, file = FaT.matrix)
# Points are the observed error rates for each consensus quality score.
# The black line shows the estimated error rates after convergence of the machine-learning algorithm.
# The red line shows the error rates expected under the nominal definition of the Q-score. 
sink("dada_learnerr_log.txt")
print("Learning forward read errors.")
errF <- learnErrors(filtFs, multithread=TRUE)
print("Learning reverse read errors.")
errR <- learnErrors(filtRs, multithread=TRUE)
sink()

pdf("ErrorPlots.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

# runs the core sample inference algorithm
sink("dada_inference_log.txt")
print("Running smaple inference algorithm for forward reads.")
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
print("Running smaple inference algorithm for reverse reads.")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
sink()
# displays dada object
dadaFs[[1]]

# Extensions: By default, the dada function processes each sample independently. 
# However, pooling information across samples can increase sensitivity to sequence 
# variants that may be present at very low frequencies in multiple samples. 
# The dada2 package offers two types of pooling. dada(..., pool=TRUE) performs 
# standard pooled processing, in which all samples are pooled together for 
# sample inference. dada(..., pool="pseudo") performs pseudo-pooling, in which 
# samples are processed independently after sharing information between samples, 
# approximating pooled sample inference in linear time.



# We now merge the forward and reverse reads together to obtain the full denoised sequences. 
# Merging is performed by aligning the denoised forward reads with the reverse-complement of
# the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. 
# By default, merged sequences are only output if the forward and reverse reads overlap by at 
# least 12 bases, and are identical to each other in the overlap region 
# (but these conditions can be changed via function arguments).

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Considerations for your own data: Most of your reads should successfully merge. 
# If that is not the case upstream parameters may need to be revisited: 
#  Did you trim away the overlap between your reads?


# Construct the sequence table ASVs
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# rows = sample names
# columns = sequence variants

# Considerations for your own data: Sequences that are much longer 
# or shorter than expected may be the result of non-specific priming.
# You can remove non-target-length sequences from your sequence table 
# (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). 
# This is analogous to “cutting a band” in-silico to get amplicons 
# of the targeted length.


# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# frequency
sum(seqtab.nochim)/sum(seqtab)


# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


# assigning taxa; must download training sets to train the algorithm, using Silva

taxa <- assignTaxonomy(seqtab.nochim, "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "./silva_species_assignment_v132.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# evaluating accuracy using mock community; must download Mock fasta
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")


# phyloseq 

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")