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

## "Phyloseq" OTU table
#ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
#               sample_data(samdf), 
#               tax_table(taxa))
F
## Construct phylogenetic tree
seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

ps <- phyloseq(tax_table(taxa), sample_data(samdf),
               otu_table(seqtab, taxa_are_rows = FALSE),phy_tree(fitGTR$tree))

ps.vibrio <- subset_taxa(ps, Genus = "Vibrio")
p1 = plot_tree(ps.vibrio, color = "Event", shape = "Group", ladderize="left") + coord_polar(theta="y")
p1