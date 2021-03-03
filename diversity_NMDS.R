rm(list=ls(all=TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

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
samdf <- data.frame(Event=sites$Site,Group=sites$Group,ID=sites$Sample) 
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
#  Define prevalence threshold as 1% of total samples
prevalenceThreshold = 0.01 * nsamples(ps)

# Filter out prevalence
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps1 = prune_taxa(keepTaxa, ps)
ps2 = tax_glom(ps1, "Genus", NArm = TRUE) #glom the pruned taxa 

top50 <- names(sort(taxa_sums(ps2), decreasing=TRUE))[1:50]
ps2.top50 <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU))
ps2.top50 <- prune_taxa(top50, ps2.top50)

# Plots from refined databases
plot_bar(ps2.top50, fill="Genus")

ordu <- ordinate(ps2.top50, method = "PCoA", distance ="bray")
p = plot_ordination(ps2.top50, ordu, color = "Event", shape = "Group")
p = p + geom_point(size=7, alpha=0.75)
p = p + scale_colour_brewer(type="qual", palette="Set1")
#p = p + geom_text(mapping = aes(label = samdf$ID), size = 4, vjust = 1.5) 

ps.vibrio <- subset_taxa(ps, Genus = "Vibrio")
