rm(list=ls(all=TRUE))

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

seqtab <- readRDS("E:/tbenhor/libraries/Documents/Oyster-Microbiota-Spring-Mortality/output/seqtab_final.rds")
taxa <- readRDS("E:/tbenhor/libraries/Documents/Oyster-Microbiota-Spring-Mortality/output/tax_final.rds")

samples.out <- rownames(seqtab)

sites <- read.csv("sites.csv", fill = FALSE, header = TRUE) 
samdf <- data.frame(Site=sites)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

#plot_richness(ps, x="Site", measures=c("Shannon", "Simpson"))

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="PCoA", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, type="" , color="Site")



