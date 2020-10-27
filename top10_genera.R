rm(list=ls(all=TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("dada2", version = "3.11")

library(dada2); packageVersion("dada2")
library(ape)
library(gridExtra)
library(knitr)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

seqtab <- readRDS("output/seqtab_cut_final.rds")
taxa <- readRDS("output/tax_cut_final.rds")

samples.out <- rownames(seqtab)

sites <- read.csv("sites_cut.csv", fill = FALSE, header = TRUE) 
samdf <- data.frame(Site=sites$Site,Group=sites$Group,ID=sites$Sample) 
rownames(samdf) <- samples.out

## "Phyloseq" object from OTU table
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

top10 <- names(sort(taxa_sums(ps2), decreasing=TRUE))[1:10]
ps2.top10 <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU))
ps2.top10 <- prune_taxa(top10, ps2.top10)

ps2.top10.BS <- subset_samples(ps2.top10, Site=="BS")
ps2.top10.SC <- subset_samples(ps2.top10, Site=="SC")
ps2.top10.NR <- subset_samples(ps2.top10, Site=="NR")

names <- taxa_names(ps2.top10.BS)

p1 = plot_heatmap(ps2.top10.BS, taxa.label = "Genus", sample.label = "Group", low="#FFFFCC", 
                  high="#000033", na.value = "white", sample.order = "Group",taxa.order = taxa_names(ps2.top10.BS),
                  trans = identity_trans())
p2 = plot_heatmap(ps2.top10.SC, taxa.label = "Genus", sample.label = "Group", low="#FFFFCC", 
                  high="#000033", na.value = "white", sample.order = "Group",taxa.order = taxa_names(ps2.top10.BS),
                  trans = identity_trans())
p3 = plot_heatmap(ps2.top10.NR, taxa.label = "Genus", sample.label = "Group", low="#FFFFCC", 
                  high="#000033", na.value = "white", sample.order = "Group",taxa.order = taxa_names(ps2.top10.BS),
                  trans = identity_trans())
grid.arrange(p1, p2, p3, nrow = 1)
