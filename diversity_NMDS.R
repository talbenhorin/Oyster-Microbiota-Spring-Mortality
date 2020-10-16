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
plot_bar(ps2.top30, fill="Genus")

ordu <- ordinate(ps2.top30, method = "PCoA", distance ="bray")
p = plot_ordination(ps2.top30, ordu, color = "Event", shape = "Group")
p = p + geom_point(size=7, alpha=0.75)
p = p + scale_colour_brewer(type="qual", palette="Set1")
p
#p = p + geom_text(mapping = aes(label = samdf$ID), size = 4, vjust = 1.5) 

ps.vibrio <- subset_taxa(ps, Genus = "Vibrio")
vtree <- rtree(ntaxa(ps.vibrio),rooted=TRUE)
physeq.vibrio = merge_phyloseq(ps.vibrio, samdf, vtree)
p1 = plot_tree(ps.vibrio, color = "Event", shape = "Group", ladderize="left") + coord_polar(theta="y")
p1