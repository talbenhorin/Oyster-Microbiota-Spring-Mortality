rm(list=ls(all=TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("DECIPHER")
BiocManager::install("phangorn")

library(ape)
library(gridExtra)
library(knitr)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(DECIPHER)
library(phangorn)
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

ps.vibrio <- subset_taxa(ps, Genus = "Vibrio")
ps.vibrio.impact <- subset_samples(ps.vibrio, Event=="Impact")
ps.vibrio.control <- subset_samples(ps.vibrio, Event=="Control")
ps.vibrio.impact.healthy <- subset_samples(ps.vibrio.impact, Group=="Healthy")
ps.vibrio.impact.erosion <- subset_samples(ps.vibrio.impact, Group=="Gill Erosion")

## Construct phylogenetic tree
seqs.control <- getSequences(ps.vibrio.control@otu_table)
seqs.impact.erosion <- getSequences(ps.vibrio.impact.erosion@otu_table)
names(seqs.control) <- seqs.control # This propagates to the tip labels of the tree
names(seqs.impact.erosion) <- seqs.impact.erosion
alignment.control <- AlignSeqs(DNAStringSet(seqs.control), anchor=NA)
alignment.impact.erosion <- AlignSeqs(DNAStringSet(seqs.impact.erosion), anchor=NA)
phang.align.control <- phyDat(as(alignment.control, "matrix"), type="DNA")
phang.align.impact.erosion <- phyDat(as(alignment.impact.erosion, "matrix"), type="DNA")
dm.control <- dist.ml(phang.align.control)
dm.impact.erosion <- dist.ml(phang.align.impact.erosion)
treeNJ.control <- NJ(dm.control) 
treeNJ.impact.erosion <- NJ(dm.impact.erosion) 
#fit = pml(treeNJ, data=phang.align)
#fitGTR <- update(fit, k=4, inv=0.2)
#fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                    rearrangement = "stochastic", control = pml.control(trace = 0))
ps.vibrio.control = merge_phyloseq(ps.vibrio.control, sample_data(samdf), treeNJ)
ps.vibrio.impact.erosion = merge_phyloseq(ps.vibrio.impact.erosion, sample_data(samdf), treeNJ.impact.erosion)

p1 = plot_tree(ps.vibrio.control, method = "treeonly", ladderize="left") #+ coord_polar(theta="y")
p2 = plot_tree(ps.vibrio.impact.erosion, method = "treeonly", ladderize="left")
grid.arrange(nrow = 1, p1, p2)

detach("package:phangorn", unload=TRUE)