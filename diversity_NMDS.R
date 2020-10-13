rm(list=ls(all=TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

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
#ps2ra = transform_sample_counts(ps2, function(x){x / sum(x)}) # Relative abundance

top30 <- names(sort(taxa_sums(ps2), decreasing=TRUE))[1:30]
ps2.top30 <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU))
ps2.top30 <- prune_taxa(top30, ps2.top30)

ordu <- ordinate(ps2.top30, method = "PCoA", distance ="bray")
p = plot_ordination(ps2.top30, ordu, color = "Event", shape = "Group")
p = p + geom_point(size=7, alpha=0.75)
p = p + scale_colour_brewer(type="qual", palette="Set1")
#p = p + geom_text(mapping = aes(label = samdf$ID), size = 4, vjust = 1.5) 

ps2.vibrio = subset_taxa(ps2.top30, Genus=="Vibrio")
#plot_bar(ps2.vibrio)
phy_tree(ps2.vibrio)
>>>>>>> b6b76ffe6a22350d170b9aac7c2dbb36cd6e147b
