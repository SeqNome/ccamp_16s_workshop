library(dada2)

setwd("16s_workshop/")

setwd("..")  
getwd()  

# Selected samples
selected_samples <- c("B1", "B3", "BW1", "BW2", "R7", "R8", "R9", "R10")

# Create subset directory
dir.create("subset_8samples", showWarnings = FALSE)

# Copy FASTQ files
for (sample in selected_samples) {
  file1 <- paste0(sample, "_R1_trim.fq.gz")
  file2 <- paste0(sample, "_R2_trim.fq.gz")
  
  if (file.exists(file1)) file.copy(file1, "subset_8samples/")
  if (file.exists(file2)) file.copy(file2, "subset_8samples/")
}

setwd("16s_workshop/data/test_out/subset_8samples/")


list.files()


r1_files <- list.files(pattern = "_R1_trim\\.fq\\.gz$")
samples <- gsub("_R1_trim\\.fq\\.gz$", "", r1_files)


# Write to samples file
writeLines(samples, "samples")

# Verify
cat("Created samples file with", length(samples), "samples:\n")
print(samples)


forward_reads <- paste0(samples, "_R1_trim.fq.gz")
reverse_reads <- paste0(samples, "_R2_trim.fq.gz")

# Create output folder for filtered reads
dir.create("filtered", showWarnings = FALSE)


# DADA2 will create filtered versions
filtered_forward_reads <- paste0("filtered/", samples, "_R1_filt.fq.gz")
filtered_reverse_reads <- paste0("filtered/", samples, "_R2_filt.fq.gz")


cat("\nChecking file availability:\n")
cat("Forward files found:", sum(file.exists(forward_reads)), "/", length(forward_reads), "\n")
cat("Reverse files found:", sum(file.exists(reverse_reads)), "/", length(reverse_reads), "\n")


# 5. QUALITY CHECK (on trimmed files)
cat("\n=== STEP 1: Quality Check ===\n")
pdf("quality_profiles_trimmed.pdf", width = 10, height = 6)
plotQualityProfile(forward_reads[1:2], aggregate = TRUE)
plotQualityProfile(reverse_reads[1:2], aggregate = TRUE)
dev.off()
cat("Quality plots saved to: quality_profiles_trimmed.pdf\n")

filtered_out <- filterAndTrim(
  fwd = forward_reads, filt = filtered_forward_reads,
  rev = reverse_reads, filt.rev = filtered_reverse_reads,
  truncLen = c(250, 200),  # Adjust based on your quality 
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE,
  verbose = TRUE
)

err_forward <- learnErrors(filtered_forward_reads, multithread = FALSE)
err_reverse <- learnErrors(filtered_reverse_reads, multithread = FALSE)


plotErrors(err_forward, nominalQ = TRUE)
plotErrors(err_reverse, nominalQ = TRUE)


derep_forward <- derepFastq(filtered_forward_reads, verbose = TRUE)
derep_reverse <- derepFastq(filtered_reverse_reads, verbose = TRUE)
names(derep_forward) <- samples
names(derep_reverse) <- samples

## ASV
dada_forward <- dada(derep_forward, err = err_forward, multithread = FALSE)
dada_reverse <- dada(derep_reverse, err = err_reverse, multithread = FALSE)



merged_amplicons <- mergePairs(
  dada_forward, derep_forward,
  dada_reverse, derep_reverse,
  verbose = TRUE
)

merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, trimOverhang=TRUE, minOverlap=170)

class(merged_amplicons) # list
length(merged_amplicons) # 
names(merged_amplicons)

class(merged_amplicons$B1) #

names(merged_amplicons$B1)

seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) 

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)



sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

summary_tab
write.table(summary_tab, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)

load("../../../data/dada2_amplicon/tax-info.RData")

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
# asv_tax <- t(sapply(tax_info, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa
# }))
# colnames(asv_tax) <- ranks
# rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)


ranks <- c("domain", "phylum", "class", "order", "family", "genus")  # No "species"

asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

# Check dimensions match
if (nrow(asv_tax) == length(asv_headers)) {
  rownames(asv_tax) <- gsub(">", "", asv_headers)
  colnames(asv_tax) <- ranks
  cat("Success! Created taxonomy table with", nrow(asv_tax), "ASVs\n")
} else {
  cat("Adjusting... Using first", min(nrow(asv_tax), length(asv_headers)), "ASVs\n")
  n <- min(nrow(asv_tax), length(asv_headers))
  asv_tax <- asv_tax[1:n, , drop = FALSE]
  rownames(asv_tax) <- gsub(">", "", asv_headers[1:n])
  colnames(asv_tax) <- ranks
}
write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

library(decontam)
packageVersion("decontam")
colnames(asv_tab)
vector_for_decontam <- c(rep(TRUE, 4), rep(FALSE, 16))
contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)
table(contam_df$contaminant)
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv", sep="\t", quote=F, col.names=NA)

library(tidyverse)
library(phyloseq)
library(vegan)
library(DESeq2)
library(dendextend)
library(viridis)

rm(list=ls())

# count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1,
#                         check.names=F, sep="\t")[ , -c(1:4)]
 #tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T,
  #                               row.names=1, check.names=F, sep="\t"))
# sample_info_tab <- read.table("sample_info.tsv", header=T, row.names=1,
#                               check.names=F, sep="\t")
# sample_info_tab$color <- as.character(sample_info_tab$color)
# sample_info_tab

library(tidyverse)
library(phyloseq)
library(vegan)
library(DESeq2)
library(dendextend)
library(viridis)


tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1, sep="\t")

sample_info_tab <- read.table("sample_info.tsv", header=T, row.names=1,
                              check.names=F, sep="\t")
# test_sample_info <- data.frame(
#   row.names = c("B1", "B2"),
#   type = c("rock", "rock"),  # Adjust based on your samples
#   char = c("altered", "glassy"),  # Or whatever characteristics
#   color = c("#FF0000", "#0000FF")  # Colors for plotting
# )

deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~type)
#deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = test_sample_info, design = ~char)

common_samples <- intersect(colnames(count_tab), rownames(sample_info_tab))
common_samples

deseq_counts <- DESeqDataSetFromMatrix(
  countData = count_tab[, common_samples],  # Subset count 
  colData = sample_info_tab[common_samples, ],  # Subset sample info
  design = ~type
)


# First estimate size factors with poscounts method
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")

# Then apply VST
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

#deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
plot(euc_clust)
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_tab$color[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols
plot(euc_dend, ylab="VST Euc. dist.")

vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues

plot_ordination(vst_physeq, vst_pcoa, color="char") + 
  geom_point(size=1) + labs(col="type") + 
  geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + 
  theme_bw() + theme(legend.position="none")

rarecurve(t(count_tab), step=100, col=sample_info_tab$color, lwd=2, ylab="ASVs", label=F)
abline(v=(min(rowSums(t(count_tab)))))

count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

plot_richness(ASV_physeq, color="char", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) +
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot_richness(ASV_physeq, x="type", color="char", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) +
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="phylum"))
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="phylum"))[,2])
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ]
class_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="class"))
class_tax_phy_tab <- tax_table(tax_glom(ASV_physeq, taxrank="class"))
phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("phylum"=phy_tmp_vec, "class"=class_tmp_vec, row.names = rows_tmp)
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$phylum == "Proteobacteria", "class"])
rownames(class_counts_tab) <- as.vector(class_tax_tab$class)
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ]
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)
identical(colSums(major_taxa_counts_tab), colSums(count_tab))
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
dim(major_taxa_proportions_tab)
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 5, ])
dim(temp_filt_major_taxa_proportions_tab)
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab
filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot)
filt_major_taxa_proportions_tab_for_plot.g <- pivot_longer(filt_major_taxa_proportions_tab_for_plot, !Major_Taxa, names_to = "Sample", values_to = "Proportion") %>% data.frame()
sample_info_for_merge<-data.frame("Sample"=row.names(sample_info_tab), "char"=sample_info_tab$char, "color"=sample_info_tab$color, stringsAsFactors=F)
filt_major_taxa_proportions_tab_for_plot.g2 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge)

ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")

ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Major_Taxa, Proportion)) +
  geom_jitter(aes(color=factor(char), shape=factor(char)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_tab_for_plot.g2$color[order(filt_major_taxa_proportions_tab_for_plot.g2$char)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank()) +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="All samples")

bw_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$type == "water"]
rock_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$type == "rock"]
filt_major_taxa_proportions_rocks_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% rock_sample_IDs, ]
filt_major_taxa_proportions_water_samples_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% bw_sample_IDs, ]

ggplot(filt_major_taxa_proportions_rocks_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) +
  geom_jitter(aes(color=factor(char), shape=factor(char)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_rocks_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_rocks_only_tab_for_plot.g$char)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Rock samples only")

ggplot(filt_major_taxa_proportions_water_samples_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) +
  geom_jitter(aes(color=factor(char)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_water_samples_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_water_samples_only_tab_for_plot.g$char)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Bottom-water samples only")

rock_sample_major_taxa_proportion_tab <- filt_major_taxa_proportions_rocks_only_tab_for_plot.g[, c(1:3)] %>% pivot_wider(names_from = Major_Taxa, values_from = Proportion) %>% column_to_rownames("Sample") %>% t() %>% data.frame()
water_sample_major_taxa_proportion_tab <- filt_major_taxa_proportions_water_samples_only_tab_for_plot.g[, c(1:3)] %>% pivot_wider(names_from = Major_Taxa, values_from = Proportion) %>% column_to_rownames("Sample") %>% t() %>% data.frame()
rock_sample_summed_major_taxa_proportions_vec <- rowSums(rock_sample_major_taxa_proportion_tab)
water_sample_summed_major_taxa_proportions_vec <- rowSums(water_sample_major_taxa_proportion_tab)
rock_sample_major_taxa_summary_tab <- data.frame("Major_Taxa"=names(rock_sample_summed_major_taxa_proportions_vec), "Proportion"=rock_sample_summed_major_taxa_proportions_vec, row.names=NULL)
water_sample_major_taxa_summary_tab <- data.frame("Major_Taxa"=names(water_sample_summed_major_taxa_proportions_vec), "Proportion"=water_sample_summed_major_taxa_proportions_vec, row.names=NULL)

ggplot(data.frame(rock_sample_major_taxa_summary_tab), aes(x="Rock samples", y=Proportion, fill=Major_Taxa)) + 
  geom_bar(width=1, stat="identity") +
  coord_polar("y") +
  scale_fill_viridis(discrete=TRUE) +
  ggtitle("Rock samples only") +
  theme_void() +
  theme(plot.title = element_text(hjust=0.5), legend.title=element_blank())

ggplot(data.frame(water_sample_major_taxa_summary_tab), aes(x="Bottom water samples", y=Proportion, fill=Major_Taxa)) + 
  geom_bar(width=1, stat="identity") +
  coord_polar("y") +
  scale_fill_viridis(discrete=TRUE) +
  ggtitle("Water samples only") +
  theme_void() +
  theme(plot.title = element_text(hjust=0.5), legend.title=element_blank())

anova(betadisper(euc_dist, sample_info_tab$type))

basalt_sample_IDs <- rock_sample_IDs[!rock_sample_IDs %in% "R7"]
basalt_euc_dist <- dist(t(vst_trans_count_tab[ , colnames(vst_trans_count_tab) %in% basalt_sample_IDs]))
basalt_sample_info_tab <- sample_info_tab[row.names(sample_info_tab) %in% basalt_sample_IDs, ]
anova(betadisper(basalt_euc_dist, basalt_sample_info_tab$char))
adonis(basalt_euc_dist~basalt_sample_info_tab$char)

basalt_vst_count_phy <- otu_table(vst_trans_count_tab[, colnames(vst_trans_count_tab) %in% basalt_sample_IDs], taxa_are_rows=T)
basalt_sample_info_tab_phy <- sample_data(basalt_sample_info_tab)
basalt_vst_physeq <- phyloseq(basalt_vst_count_phy, basalt_sample_info_tab_phy)
basalt_vst_pcoa <- ordinate(basalt_vst_physeq, method="MDS", distance="euclidean")
basalt_eigen_vals <- basalt_vst_pcoa$values$Eigenvalues

plot_ordination(basalt_vst_physeq, basalt_vst_pcoa, color="char") + 
  labs(col="type") + geom_point(size=1) + 
  geom_text(aes(label=rownames(basalt_sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=25, y=68, label="Highly altered vs glassy") +
  annotate("text", x=25, y=62, label="Permutational ANOVA = 0.003") + 
  coord_fixed(sqrt(basalt_eigen_vals[2]/basalt_eigen_vals[1])) + ggtitle("PCoA - basalts only") + 
  scale_color_manual(values=unique(basalt_sample_info_tab$color[order(basalt_sample_info_tab$char)])) + 
  theme_bw() + theme(legend.position="none")

basalt_count_phy <- otu_table(count_tab[, colnames(count_tab) %in% basalt_sample_IDs], taxa_are_rows=T)
basalt_count_physeq <- phyloseq(basalt_count_phy, basalt_sample_info_tab_phy)
basalt_deseq <- phyloseq_to_deseq2(basalt_count_physeq, ~char)
basalt_deseq <- DESeq(basalt_deseq)

deseq_res_altered_vs_glassy <- results(basalt_deseq, alpha=0.01, contrast=c("char", "altered", "glassy"))
summary(deseq_res_altered_vs_glassy)
sigtab_res_deseq_altered_vs_glassy <- deseq_res_altered_vs_glassy[which(deseq_res_altered_vs_glassy$padj < 0.01), ]
summary(sigtab_res_deseq_altered_vs_glassy)
sigtab_deseq_altered_vs_glassy_with_tax <- cbind(as(sigtab_res_deseq_altered_vs_glassy, "data.frame"), as(tax_table(ASV_physeq)[row.names(sigtab_res_deseq_altered_vs_glassy), ], "matrix"))
sigtab_deseq_altered_vs_glassy_with_tax[order(sigtab_deseq_altered_vs_glassy_with_tax$baseMean, decreasing=T), ]