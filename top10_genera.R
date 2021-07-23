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
library(scales)
theme_set(theme_bw())

seqtab <- readRDS("output/seqtab.rds")
taxa <- readRDS("output/taxtab.rds")

samples.out <- rownames(seqtab)

sites <- read.csv("sites_all.csv", fill = FALSE, header = TRUE) 
samdf <- data.frame(Stage=sites$Stage,ID=sites$Sample,OV=sites$OV) 
rownames(samdf) <- samples.out

## "Phyloseq" object from OTU table
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

ps <- subset_samples(ps, OV=="Yes")

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

top20 <- names(sort(taxa_sums(ps2), decreasing=TRUE))[1:20]
ps2.top20 <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU))
ps2.top20 <- prune_taxa(top20, ps2.top20)

# plot_bar(ps2.top20, fill="Genus")

ordu <- ordinate(ps2.top20,"NMDS","bray")
p = plot_ordination(ps2.top20, ordu, color = "Stage")
p = p + geom_point(size=4, alpha=0.75)
p = p + scale_colour_brewer(type="qual", palette="Set1")

ps.A <- subset_samples(ps2.top20, Stage=="Female")
ps.B <- subset_samples(ps2.top20, Stage=="Oligo Female")
ps.C <- subset_samples(ps2.top20, Stage=="Virilescent Female")

p1 = plot_heatmap(ps.A, taxa.label = "Genus", sample.label = "ID",low="white", high="#000033", 
                  na.value = "white", sample.order = "ID",taxa.order = taxa_names(ps.A),
                  trans = identity_trans())
p1 = p1 + theme(axis.text.x = element_text(size=7, angle=0, hjust=0.5, vjust=0.95)) +
                  theme(axis.title.y=element_blank()) +                
                  theme(legend.position="none") +
                  ggtitle("Female") +
                  theme(axis.title.x=element_blank())
p2 = plot_heatmap(ps.B, taxa.label = "Genus", sample.label = "ID",low="white", high="#000033", 
                  na.value = "white", sample.order = "ID",taxa.order = taxa_names(ps.B),
                  trans = identity_trans())
p2 = p2 + theme(axis.text.x = element_text(size=7, angle=0, hjust=0.5, vjust=0.95)) +
                  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) + 
                  theme(legend.position="none") +
                  ggtitle("Oligo Female") +
                  theme(axis.title.x=element_blank())
p3 = plot_heatmap(ps.C, taxa.label = "Genus", sample.label = "ID",low="white", high="#000033", 
                  na.value = "white", sample.order = "ID",taxa.order = taxa_names(ps.C),
                  trans = identity_trans())
p3 = p3 + theme(axis.text.x = element_text(size=7, angle=0, hjust=0.5, vjust=0.95),
                  legend.title = element_text(size = 16)) +
                  labs(fill = "Relative\nabundance") +
                  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +
                  ggtitle("Virilescent Female") +
                  theme(axis.title.x=element_blank())
grid.arrange(p1, p2, p3, nrow = 1)
