#Expression analysis for trancriptomic data
#Aung 16Feb 2026
##############################################################
# 1) LOAD REQUIRED PACKAGES
##############################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") # allows install of packages from Bioconducter project 
BiocManager::install(c("affy", "limma", "hgu133a.db", "annotate"), version = '3.22')

library(affy) # methods for Affymetrix Oligonucleotide Arrays
library(limma) # Linear Models for Microarray and Omics Data
library(hgu133a.db) # Affymetrix Affymetrix HG-U133A Array annotation data (chip hgu133a)
library(annotate) # Annotate for microarrays
library(ggplot2) # Grammar of graphics, tidyverse, creates good visualizations
library(pheatmap) # Pretty Heatmaps


##############################################################
# 2) READ RAW CEL FILES + RMA NORMALIZATION
##############################################################
untar("/BI428-Transciptomics-Analysis/GSE5389_RAW.tar", # Extract files from or list the contents of a tar archive
      exdir = "/BI428-Transciptomics-Analysis/GSE5389_RAW/") # The directory to extract files to

# A TAR file (Tape Archive file) is an archive format used to package multiple files 
# together for easier storage and sharing. 
# Unlike ZIP files, TAR files are not compressed by default

# folder that contains ONLY the CEL files
cel_path <- "/BI428-Transciptomics-Analysis/GSE5389_RAW/"

raw <- ReadAffy(celfile.path = cel_path) # Read CEL files into an Affybatch
# celfile.path = a character denoting the path ReadAffy should look for cel files

# What is an Affybatch?
# This is a class representation for Affymetrix GeneChip probe level data. 
# The main component are the intensities from multiple arrays of the same CDF type

# Perform RMA (log2, background-corrected, quantile-normalized)
eset_rma <- rma(raw) # Robust Multi-Array Average expression measure

# This function converts an AffyBatch object into an ExpressionSet 
# object using the robust multi-array average (RMA) expression measure.

# What is an ExpressionSet object?
# Container for high-throughput assays and experimental metadata. ExpressionSet class is derived from eSet, 
# and requires a matrix named exprs as assayData member.

expr_rma <- exprs(eset_rma)   # matrix (probe x sample)

# exprs = hese generic functions access the expression and error measurements 
# of assay data stored in an object derived from the eSet-class

##############################################################
# 3) EXTRACT SAMPLE NAMES + CREATE BD/CTL LABELS
##############################################################

samples <- colnames(expr_rma) 

# Create Dx vector manually using GSM identifiers
# Modify based on your sample list:
# Example (you MUST verify!)
bd_ids  <- c("GSM123243.cel.gz","GSM123244.cel.gz","GSM123245.cel.gz","GSM123246.cel.gz","GSM123247.cel.gz",
             "GSM123248.cel.gz","GSM123249.cel.gz","GSM123250.cel.gz","GSM123251.cel.gz","GSM123252.cel.gz")

ctl_ids <- c("GSM123253.cel.gz","GSM123254.cel.gz","GSM123255.cel.gz","GSM123256.cel.gz","GSM123257.cel.gz",
             "GSM123258.cel.gz","GSM123259.cel.gz","GSM123260.cel.gz","GSM123261.cel.gz","GSM123262.cel.gz","GSM123263.cel.gz")

# c() means combine, used for creating vectors 

Dx <- ifelse(samples %in% bd_ids, "BD", 
             ifelse(samples %in% ctl_ids, "CTL", NA))

Dx <- factor(Dx, levels = c("CTL","BD")) # changes Dx to a factor

##############################################################
# 4) PROBE → GENE SYMBOL MAPPING + COLLAPSE PROBES
##############################################################

probe_ids <- rownames(expr_rma)

symbols <- getSYMBOL(probe_ids, "hgu133a") # getSYMBOL(x, data)

expr2 <- expr_rma
rownames(expr2) <- symbols

# Remove NAs
expr2 <- expr2[!is.na(rownames(expr2)), ]

# Collapse probes (median)
expr_gene <- avereps(expr2, ID = rownames(expr2), FUN = median)
# Average Over Irregular Replicate Probes
# Condense a microarray data object so that values 
# for within-array replicate probes are replaced with their average.
# ID = probe identifier

##############################################################
# 5) LIMMA DIFFERENTIAL EXPRESSION (BD vs CTL)
##############################################################

design <- model.matrix(~ Dx)
# model.matrix creates a design (or model) matrix, e.g., by expanding factors 
# to a set of dummy variables (depending on the contrasts) and expanding interactions similarly
fit <- lmFit(expr_gene, design) # Fit linear model for each gene given a series of arrays
fit <- eBayes(fit)
# Given a linear model fit from lmFit, compute moderated t-statistics, moderated 
# F-statistic, and log-odds of differential expression by empirical Bayes moderation 
# of the standard errors towards a global value

de_all <- topTable(fit, coef = "DxBD", number = Inf)
# Extract a table of the top-ranked genes from a linear model fit
# coef = column number or column name specifying which coefficient or contrast of the linear model is of interest
# maximum number of genes to list

# logFC = estimate of the log2-fold-change corresponding to the effect or contrast
# AveExpr = average log2-expression for the probe over all arrays and channels
# t = moderated t-statistic
# P.Value = raw p-value
# adj.P.Val = adjusted p-value or q-value
# B = log-odds that the gene is differentially expressed

##############################################################
# SAVE RESULTS
##############################################################

write.csv(de_all, "GSE5389_DE_results_RMA.csv", row.names = TRUE)

library(readr)
GSE5389_DE_results_RMA <- read_csv("GSE5389_DE_results_RMA.csv")
View(GSE5389_DE_results_RMA)

##############################################################
# PCA PLOT (BD vs CTL)
##############################################################

# Transpose so samples = rows, genes = columns
pca <- prcomp(t(expr_gene), scale. = TRUE)
# Performs a principal components analysis on the given data 
# matrix and returns the results as an object of class prcomp
# t() is returing the transpose of expr_gene
# scale. = a logical value indicating whether the variables 
# should be scaled to have unit variance before the analysis takes place


# Create PCA dataframe
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Dx = Dx
)

# Variance explained
var_explained <- round(100 * summary(pca)$importance[2, 1:2], 2)
# rounding the numbers to 2 decimal places

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Dx)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA Plot (RMA normalized expression)",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  scale_color_manual(values = c("CTL" = "blue", "BD" = "red"))

##############################################################
# VOLCANO PLOT
##############################################################

volc <- de_all
volc$Gene <- rownames(volc)
volc$negLogP <- -log10(volc$P.Value) # -log10 of the P.Value

logFC_cutoff <- 1 # why is the cutoff 1??
p_cutoff <- 0.05

volc$Sig <- with(volc, # Evaluate an R expression in an environment constructed from data, possibly modifying (a copy of) the original data.
                 ifelse(P.Value < p_cutoff & logFC >  logFC_cutoff, "Up", 
                        # is the P.Value is less than (<) 0.05 & if the logFC is greater than (>) 1, the gene is unregulated
                        ifelse(P.Value < p_cutoff & logFC < -logFC_cutoff, "Down", "NS"))
                        # is the P.Value is less than (<) 0.05 & is the logFC is less than (<) -1, the gene is downregulated
                        # if neither of these conditions are true, the gene is not significant
)
library(tidyverse)

volc |>
  filter(Sig == "Down")

# "Up" means unregulated, "Down" means down regulated, "NS" means not significant

ggplot(volc, aes(logFC, negLogP, color = Sig)) +
  geom_point(alpha = 0.7, size = 1.8) + # alpha = opacity - range from 0 to 1, size = size of points and text 
  scale_color_manual(values=c("Up"="red","Down"="blue","NS"="darkgreen")) +
  geom_vline(xintercept=c(-logFC_cutoff, logFC_cutoff),
             linetype="dashed", col="black") +
  geom_hline(yintercept= -log10(p_cutoff),
             linetype="dashed", col="black") +
  theme_minimal(base_size=14) + # base font size, given in pts
  labs(title="Volcano Plot: BD vs Control (GSE5389, RMA-normalized)",
       x="log2 Fold Change",
       y="-log10(P-value)")

ggsave("Volcano_GSE5389_RMA.png", width=8, height=5)


##############################################################
# HEATMAP OF TOP 50 DE GENES
##############################################################

# Select top 50 DE genes
top_genes <- rownames(de_all)[1:50] # taking the top 50 genes out of de_all

# Subset expression matrix
expr_top <- expr_gene[top_genes, ]

# Z-score scale per gene
expr_scaled <- t(scale(t(expr_top)))

# Annotation for columns
annotation_col <- data.frame(Diagnosis = Dx)
rownames(annotation_col) <- colnames(expr_scaled)

# Heatmap
pheatmap(expr_scaled,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 6, # I changed the font size to make the plot more readable
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 50 Differentially Expressed Genes (RMA)")


ggsave("PCA_GSE5389.png", width = 7, height = 6, dpi = 300)
png("Heatmap_GSE5389.png", width = 1200, height = 1400, res = 150)
pheatmap(expr_scaled,
         annotation_col = annotation_col,
         main = "Top 50 DE Genes")
dev.off()


##############################################################
# HEATMAP FOR SPECIFIC PATHWAY GENES
##############################################################

pathway <- c("KEAP1","NFE2L2","HMOX1","NQO1","TXNRD1",
             "GCLC","GCLM","GSR","SRXN1","SLC7A11",
             "GPX2","PRDX1","GSTA1","GSTM1","GSTM2","GSTP1")

# Filter only genes that actually exist in expression matrix
genes_present <- pathway[pathway %in% rownames(expr_gene)]
genes_missing <- pathway[!pathway %in% rownames(expr_gene)]

# Outputs the objects, concatenating the representations 
# cat performs much less conversion than print
cat("Genes found in matrix:\n")
print(genes_present)

cat("Genes NOT found in matrix:\n")
print(genes_missing)

# Subset expression matrix
expr_subset <- expr_gene[genes_present, ]

# Z‑score transform by gene
expr_scaled <- t(scale(t(expr_subset)))

# Annotation for samples
annotation_col <- data.frame(Diagnosis = Dx)
rownames(annotation_col) <- colnames(expr_scaled)

# Heatmap
pheatmap(expr_scaled,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 10,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "NRF2 / KEAP1 Pathway Gene Expression (RMA)")
``


png("NRF2_Pathway_Heatmap.png", width = 1100, height = 1100, res = 150)
pheatmap(expr_scaled,
         annotation_col = annotation_col,
         main = "NRF2 / KEAP1 Pathway Gene Expression")
dev.off()

