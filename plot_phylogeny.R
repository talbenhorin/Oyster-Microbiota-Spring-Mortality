rm(list=ls(all=TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("DECIPHER")
BiocManager::install("phangorn")
BiocManager::install("dada2", version = "3.11")

library(dada2); packageVersion("dada2")
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
samdf <- data.frame(Site=sites$Site,ID=sites$Sample,Event=sites$Mort,Group=sites$Group) 
rownames(samdf) <- samples.out

## "Phyloseq" object from OTU table
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

ps.vibrio = prune_taxa(keepTaxa, ps)

## Phylogenetic Trees by Genus
ps.vibrio <- subset_taxa(ps, Genus=="Vibrio")
## Construct phylogenetic tree
seqs <- getSequences(ps.vibrio@otu_table)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
ps.vibrio = merge_phyloseq(ps.vibrio, sample_data(samdf), treeNJ)

p1 = plot_tree(ps.vibrio, color = "Site", shape = "Group", ladderize = TRUE) 

detach("package:phangorn", unload=TRUE)