rm(list=ls(all=TRUE))

library(gridExtra)
library(knitr)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

seqtab <- readRDS("E:/tbenhor/libraries/Documents/Oyster-Microbiota-Spring-Mortality/output/seqtab_water_final.rds")
taxa <- readRDS("E:/tbenhor/libraries/Documents/Oyster-Microbiota-Spring-Mortality/output/tax_water_final.rds")

samples.out <- rownames(seqtab)

sites <- read.csv("sites_water.csv", fill = FALSE, header = TRUE) 
samdf <- data.frame(Site=sites$Site,Event=sites$Event) #Fake data b/c the plot_ordination plot-by bug
rownames(samdf) <- samples.out

## "Phyloseq" OTU table
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

## Filtering
#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)

# Filter out prevalence
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps1 = prune_taxa(keepTaxa, ps)
ps2 = tax_glom(ps1, "Genus", NArm = TRUE) #glom the pruned taxa 

## Abundance value transformation#
ps2ra = transform_sample_counts(ps2, function(x){x / sum(x)}) # Relative abundance

top20 <- names(sort(taxa_sums(ps2ra), decreasing=TRUE))[1:20]
ps2ra.top20 <- transform_sample_counts(ps2ra, function(OTU) OTU/sum(OTU))
ps2ra.top20 <- prune_taxa(top20, ps2ra.top20)
plot_bar(ps2ra.top20, fill="Genus")

## Ordination Plots from transformed data
ordu <- ordinate(ps2ra.top20, method = "PCoA", distance ="bray")
plot_ordination(ps2ra.top20, ordu, color = "Event")

