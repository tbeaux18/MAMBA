#!/usr/bin/env Rscript
library(dada2)
library(tidyverse)
library(magrittr)

# declaring path to trimmed fastq files
fastq_path <- '/homes/FALL2019SHARED/TrimmedProjectData'


# converting to datafraem to merge into experimental design file 
trimmed_df <- as.data.frame(list.files(fastq_path)) 
colnames(trimmed_df) <- c("fastq_path")
trimmed_df %<>% mutate(sample_name = str_extract(fastq_path, pattern = "MB\\d{3}"), read_pair = str_extract(fastq_path, pattern = 'R1|R2')) %<>% spread(read_pair, fastq_path)

# merging the experimental design file with new trimmed paths
exp_data <- read.csv('/homes/tbaker8/dada_analysis/20191015_metagenomics_exp_design.csv', stringsAsFactors = F, header = T)
exp_data %<>% inner_join(., trimmed_df, by=c("sample_name"="sample_name"))


# isolating mock and blank data to run through first
mock_data <- exp_data %>% filter(replicate_number == 'Mock' | replicate_number == 'Blank')
forward_mock <- paste(fastq_path, mock_data$R1, sep="/")
reverse_mock <- paste(fastq_path, mock_data$R2, sep="/")
mock_names <- as.vector(as.character(mock_data$sample_name))


# use this to go through all samples
# forward_reads <- sort(list.files(fastq_path, pattern="_R1_001_trimmed.fastq", full.names = TRUE))
# reverse_reads <- sort(list.files(fastq_path, pattern="_R2_001_trimmed.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

pdf("mock_quality_profile.pdf")
for(i in 1:length(forward_mock)){
  
  fp <- plotQualityProfile(forward_mock[i])
  rp <- plotQualityProfile(reverse_mock[i])
  print(fp)
  print(rp)
}
dev.off()




dada_path <- "~/dada_analysis"
# Place filtered files in filtered/ subdirectory
foward_filtered <- file.path(dada_path, "filtered", paste0(mock_names, "_F_filt.fastq.gz"))
reverse_filtered <- file.path(dada_path, "filtered", paste0(mock_names, "_R_filt.fastq.gz"))
names(foward_filtered) <- mock_names
names(reverse_filtered) <- mock_names


# figures out usuable reads
# need to determine aggregate where to truncate
out <- filterAndTrim(forward_mock, foward_filtered, reverse_mock, reverse_filtered, 
                     truncLen=c(220,200), maxN=0, maxEE=c(2,2), truncQ=2, 
                     rm.phix=TRUE, compress=TRUE, multithread=TRUE)
dout <- as.data.frame(out)
FaT.matrix <- file.path(dada_path, "mock_filtered_trimmed_reads.csv")
write.csv(dout, file = FaT.matrix)
# dout %<>% rownames_to_column("file_path")


# dereplication after trimming condenses all same reads to increase computational efficiency
# common in ASV inference
forward_derep <- derepFastq(foward_filtered, verbose=TRUE)
names(forward_derep) <- mock_names
reverse_derep <- derepFastq(reverse_filtered, verbose=TRUE)
names(reverse_derep) <- mock_names


forward_errors <- learnErrors(forward_derep, multithread=TRUE, randomize=TRUE)
reverse_errors <- learnErrors(reverse_derep, multithread=TRUE, randomize=TRUE)


# pdf("ErrorPlots.pdf")
# par(mfrow=c(1,2))
plotErrors(forward_errors, nominalQ=TRUE)
plotErrors(reverse_errors, nominalQ=TRUE)
# dev.off()


# runs the core sample inference algorithm
forward_dada <- dada(forward_derep, err=forward_errors, multithread=TRUE)
reverse_dada <- dada(reverse_derep, err=reverse_errors, multithread=TRUE)

# displays dada object
dadaFs[[1]]


# Merges read 1 and read 2 into a contig if overlapping by 12 base pairs by default
# Need to change this overlap
mergers <- mergePairs(forward_dada, forward_derep, reverse_dada, reverse_derep, minOverlap=150, trimOverhang=TRUE, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


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



seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# frequency 25%
chimera.freq <- 1-(sum(seqtab.nochim)/sum(seqtab))


getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=mock_names, dada2_input=dout[,1],
                          filtered=dout[,2], dada_f=sapply(forward_dada, getN),
                          dada_r=sapply(reverse_dada, getN), merged=sapply(mergers, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/dout[,1]*100, 1))
summary_tab %<>% rownames_to_column("sample_name")
mock_data %<>% inner_join(.,summary_tab, by=c("sample_name"="sample_name"))
# write.csv(summary_tab, file="summarytab.csv")


ggplot(mock_data, aes(x=as.factor(sample_name), y=dada2_input)) + 
  geom_bar(stat="identity", position="dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~replicate_number)


# assigning taxa; must download training sets to train the algorithm, using Silva
taxa <- assignTaxonomy(seqtab.nochim, "~/dada_analysis/silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC=T)
taxa <- addSpecies(taxa, "~/dada_analysis/silva_species_assignment_v132.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)



# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
asv_map <- data.frame(asv_headers, asv_seqs)
asv_map %<>% mutate(asv_headers = str_sub(asv_headers, 2, -1))
asv_tax_new <- as.data.frame(asv_tax) %>% rownames_to_column("asv_names")

asv_map_total <- inner_join(asv_map, asv_tax_new, by=c("asv_headers"="asv_names"))

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "mock_ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "mock_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "mock_ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)


library("decontam")
colnames(asv_tab)
vector_samples <- colnames(asv_tab)
vector_for_decontam <- c('MB001', 'MB007', 'MB013', 'MB019', 'MB025', 'MB031', 'MB037', 'MB043', 'MB049', 'MB061', 'MB067', 'MB073', 'MB079', 'MB085', 'MB091', 'MB097', 'MB103', 'MB109', 'MB115', 'MB121', 'MB127', 'MB133', 'MB139', 'MB145', 'MB151', 'MB157', 'MB163', 'MB169', 'MB175', 'MB181', 'MB187', 'MB193', 'MB199', 'MB205', 'MB211', 'MB217', 'MB223', 'MB229', 'MB235', 'MB241', 'MB247', 'MB253', 'MB259', 'MB265', 'MB271')
vector_for_decontam_idx <- vector_samples %in% vector_for_decontam
contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam_idx)

table(contam_df$contaminant) # identified 6 as contaminants

# getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

asv_tax[row.names(asv_tax) %in% contam_asvs, ]

# making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

# making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

# making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

## and now writing them out to files
write(asv_fasta_no_contam, "mock_ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "mock_ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "mock_ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)



library(Biostrings)
fasta_file = readDNAStringSet("ssrrna_ref.fasta")
seq_name = names(fasta_file)
sequence = paste(fasta_file)
fasta_df <- data.frame(seq_name, sequence)
fasta_df %<>% mutate(org_name = str_extract(seq_name, pattern = '.+?(?=_\\d)')) %<>% filter(!org_name %in% c('Cryptococcus', 'Saccharomyces'))

org_name <- unique(fasta_df$org_name)

# mock.ref <- getSequences(file.path(dada_path, "zymo_ref1.fasta"))
# mock.ref <- getSequences(file.path(dada_path, "ssrrna_ref.fasta"))
mock_vector <- c('MB002', 'MB008', 'MB014', 'MB020', 'MB026', 'MB032', 'MB038', 'MB044', 'MB050', 'MB062', 'MB068', 'MB074', 'MB080', 'MB086', 'MB092', 'MB098', 'MB104', 'MB110', 'MB116', 'MB122', 'MB128', 'MB134', 'MB140', 'MB146', 'MB152', 'MB158', 'MB164', 'MB170', 'MB176', 'MB182', 'MB188', 'MB194', 'MB200', 'MB206', 'MB212', 'MB218', 'MB224', 'MB230', 'MB236', 'MB242', 'MB248', 'MB254', 'MB260', 'MB266', 'MB272')
unqs.mock <- seqtab.nochim[mock_vector,]

mock.total <- mock_data %>% select(sample_name, dada2_input)
# mock.norm <- unqs.mock inner_join(unqs.mock, mock.total, by

sum_total <- setNames(data.frame(matrix(ncol = 2, nrow = length(mock_vector))), c("community_sum", "sum_exact_matches"))
mock_sum <- setNames(data.frame(matrix(ncol = 8, nrow = length(mock_vector))), c(org_name))
rownames(mock_sum) <- mock_vector

datalist = list()

for (i in mock_vector){
  mockf <- unqs.mock[i,]
  mock_sorted <- sort(mockf[mockf>0], decreasing=TRUE)
  sum_total[i, "community_sum"] <- length(mock_sorted)
  
  if (sum_total[i, "community_sum"] != 0) {
    sumthing <- sapply(names(mock_sorted), function(x) any(grepl(x, fasta_df$sequence)))
    sum_df <- data.frame(match = sumthing) %>% rownames_to_column("mock_seqs")
    sum_df$sample_name <- i

    join_df <- inner_join(sum_df, asv_map_total, by=c("mock_seqs"="asv_seqs"))
    datalist[[i]] <- join_df
    
    match.ref <- sum(sapply(names(mock_sorted), function(x) any(grepl(x, fasta_df$sequence))))
    
    sum_total[i, "sum_exact_matches"] <- sum(match.ref)
  } else {
    sum_total[i, "sum_exact_matches"] <- 0
  }
}

big_data = do.call(rbind, datalist)

big_data %<>% filter(match == 'TRUE')

try <- big_data %>% group_by(sample_name, Family) %>% summarise(n=n())




# phyloseq
library(ggplot2)
library(phyloseq)
library(Biostrings)

# manipulate dataframe to pass in as samdf
samdf <- mock_data %>% select(sample_name, replicate_number, analyst_number, extraction_kit)
rownames(samdf) <- mock_data$sample_name

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))


# ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


plot_richness(ps)

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")


top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="extraction_kit", fill = "Family") + facet_wrap(~replicate_number)


organisms <- c("Pseudomonas aeruginosa", "Escherichia coli", "Salmonella enterica", "Lactobacillus fermentum", "Enterococcus faecalis", "Staphylococcus aureus", "Listeria monocytogenes", "Bacillus subtilis") 
comp_pct <- c(4.2, 10.1, 10.4, 18.4, 9.9, 15.5, 14.1, 17.4)
cp_num <- c(4,7,7,5,4,6,6,10)
stdcomp <- data.frame(organisms, comp_pct, cp_num)


library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")

count_tab <- asv_tab[,mock_vector]
tax_tab <- asv_tax
mock_only_data <- mock_data[mock_data$sample_name %in% mock_vector,]

# apply(count_tab, 2, function(x) colSums(x))

# first we need to make a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = mock_only_data, design = ~extraction_kit) 
# we have to include the "colData" and "design" arguments because they are 
# required, as they are needed for further downstream processing by DESeq2, 
# but for our purposes of simply transforming the data right now, they don't 
# matter
deseq_counts <- deseq_counts[rowSums(counts(deseq_counts)) > 1,]

deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")

deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts2)
# NOTE: If you get this error here with your dataset: "Error in
# estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means", that
# can be because the count table is sparse with many zeroes, which is common
# with marker-gene surveys. In that case you'd need to use a specific
# function first that is equipped to deal with that. You could run:
# deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
# deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))



rarecurve(t(count_tab), step=100, lwd=2, ylab="ASVs", label=F)
abline(v=(min(rowSums(t(count_tab)))))
