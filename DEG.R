# -------------------------------------------#
# Section 0: Load Required Libraries---
# -------------------------------------------#
library(limma)
library(edgeR)
library(dplyr)
library(stringr)
library(pheatmap)
library(ggplot2)

# -------------------------------------------#
#Section 1: read and explore the data count----
# -------------------------------------------#
# Define the path to the featureCounts_output.txt file
file_path <- "./featureCounts_output.txt"

# Read the file into R
featureCounts_data <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1)

# View the first few rows of the data
head(featureCounts_data)

# Modify the column names by removing the unwanted part of the filenames
colnames(featureCounts_data) <- gsub("(_1)?_trimmed_Aligned.sortedByCoord.out.bam", "", colnames(featureCounts_data))

# View the first few rows to check the output
head(featureCounts_data)

# View the summary of the data
summary(featureCounts_data)

# View the structure of the data
str(featureCounts_data)

# -------------------------------------------#
#Section 2: Prepare the Metadata, gene annotation and gene symbole----
# -------------------------------------------#
## 1: Metadata Setup----
# Define the metadata (treatments, replicate) for the samples
sample_metadata <- data.frame(
  row.names = colnames(featureCounts_data),
  treatments = c("Mature biofilm", "Mature biofilm", "Mature biofilm", "Thin biofilm", "Thin biofilm", "Thin biofilm", "Early biofilm", "Early biofilm", "Early biofilm")
)
# Check the metadata 
head(sample_metadata)

# Create DGEList object
dge <- DGEList(counts = featureCounts_data, group = sample_metadata$treatments)

# Check the DGEList object
dge

# Ensure the sample names (row names in sample_metadata) match column names in featureCounts_data
colnames(featureCounts_data) <- rownames(sample_metadata)

# Assign the treatments to dge$samples$group
dge$samples$group <- sample_metadata$treatments

# Check the updated dge$samples to ensure the treatments are correctly added
head(dge$samples)

## 2: Add Gene Annotation -----
# Load the RData file into R that contains the gene annotation information
load("fc2.RData")

# Check the objects in the RData file
ls()

# View the annotation data
head(fc2$annotation)

# Extract annotation into a separate variable
gene_annotation <- fc2$annotation

# View the extracted annotation
head(gene_annotation)

# Save the dge list with gene annotation
save(dge, file = "dge_with_gene_annotation.RData")
 
## 3: Add Gene Symbols----
# Define the path to the GTF file
gtf_file_path <- "./genomic.gtf"

# Read the GTF file into R
gtf_data <- read.delim(gtf_file_path, header = FALSE, comment.char = "#", sep = "\t")

# Check the first few rows to understand its structure
head(gtf_data)

# Extract GeneID from the V9 column

# Add new columns for gene_id and gene_symbol
gtf_data <- gtf_data %>%
  mutate(
    gene_id = str_extract(V9, "gene_id\\s+\\w+") %>% str_remove("gene_id\\s+"),
    gene_symbol = str_extract(V9, "gene\\s+\\w+") %>% str_remove("gene\\s+")
  )

# Ensure that both Gene_ID and gene_id columns are correctly named and cleaned
# Clean the gene_id in gtf_data and Gene_ID in gene_annotation 
gtf_data$gene_id <- trimws(gtf_data$gene_id)
gene_annotation$GeneID <- trimws(gene_annotation$GeneID)

# Check the column names in both data frames
colnames(gene_annotation)  # Ensure Gene_ID is the correct name for gene ID in gene_annotation
colnames(gtf_data)        # Ensure gene_id is the correct name for gene ID in gtf_data

# Rename columns for consistency
colnames(gene_annotation)[colnames(gene_annotation) == "GeneID"] <- "gene_id"

# Check the column names 
colnames(gene_annotation)  

# Merge the gene_symbol from gtf_data to gene_annotation based on matching gene_id 
gene_annotation <- merge(gene_annotation, gtf_data[, c("gene_id", "gene_symbol")], 
                         by = "gene_id", all.x = TRUE)
# Check the updated gene_annotation
head(gene_annotation)

# View the extracted GeneID and GeneSymbol
head(gtf_data[, c("gene_id", "gene_symbol")])

# Add the gene annotation information to the DGEList object
dge$genes <- gene_annotation[match(rownames(dge), gene_annotation$gene_id), ]

head(dge$genes, 100)

# Replace NA values in the gene_symbol column with "unknown"
dge$genes$gene_symbol[is.na(dge$genes$gene_symbol)] <- "unknown"

# Count number of unique values of gene_id in gtf_data
length(unique(dge$genes$gene_id)) #6470
length(unique(gene_annotation$gene_id)) #6470
length(unique(dge$genes$gene_symbol)) # 5220

# Check the updated DGEList object with Gene Symbols
head(dge$genes, 100)

# -------------------------------------------#
# Section 3: Filtering low expression data----
# -------------------------------------------#
# 1- Order genes by the sum of counts across all samples (decreasing)
o <- order(rowSums(dge$counts), decreasing = TRUE)
dge <- dge[o,]

# View the first few rows after sorting
head(dge$counts)

# 2-  Remove duplicate genes (based on gene symbol)

# Check for duplicated gene symbols
d <- duplicated(dge$genes$gene_symbol)
length(dge$genes$gene_symbol) #6470
head(d)
table(d)  #count of gene symbols TRUE (duplicates = 1250) and FALSE (unique= 5220) values

# Remove the duplicated entries
dge <- dge[!d,]
length(dge$genes$gene_symbol) #5220
nrow(dge) #5220 genes

# View the updated DGEList object
head(dge$genes, 50)

# 3-filter genes with low expressions
# Add grouping information to the DGEList 
dge$samples$group <- as.factor(c("Mature biofilm", "Mature biofilm", "Mature biofilm", "Thin biofilm", "Thin biofilm", "Thin biofilm", "Early biofilm", "Early biofilm", "Early biofilm"))

# Use filterByExpr to filter genes with low expression
keep <- filterByExpr(dge) #min.count = 10


# View the number of genes kept (TRUE means kept)
table(keep)

#kept = 5078, removed = 142   

# Subset the DGEList to keep only the genes with high expression
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Check the number of genes retained after filttering
#5220 before filtering
nrow(dge$genes)#5078  after filtering4

# -------------------------------------------#
# Section 4: Normalizing Counts and exploring variation and separation of the data----
# -------------------------------------------#
# -------------------------------------------#
## Normalizing Counts----
# -------------------------------------------#
## Normalize filtered data using TMM (Trimmed Mean of M-values) method----
dge <- calcNormFactors(dge, method = "TMM")

# Display normalization factors
dge$samples$norm.factors
#0.8024760 0.7409500 0.6234132 1.9016864 1.6059576 0.7229335 1.0077664 1.1738045 1.03294235

dge$samples
#-------------------------------------#
## Multidimensional Scaling -----
#----------------------------------------#
# Create a color mapping for the treatments
treatment_colors <- c("Mature biofilm" = "blue", "Thin biofilm" = "red", "Early biofilm" = "darkgreen")

# Run the MDS plot 
mds <- plotMDS(dge, 
               col = treatment_colors[dge$samples$group],  # Color by treatment group
               pch = 19,                                  # Point shape
               cex = 1.5,                                 # Point size
               main = "MDS Plot by Treatment",             # Plot title
               plot = FALSE)                               # Don't plot immediately

# MDS plot with treatments in different colors
plotMDS(dge, col = treatment_colors[dge$samples$group], 
        main = "MDS Plot by Treatment", 
        xlab = "Leading logFC dim 1", 
        ylab = "Leading logFC dim 2")

# Add a legend
legend("bottomright", legend = names(treatment_colors), 
       fill = treatment_colors, title = "Treatment", bty = "n")

#------------------------------------------------------------#
## Estimate dispersion- variation of genes across replicates----
#--------------------------------------------------------------#
# Define the design matrix 
design <- model.matrix(~0 + sample_metadata$treatments)

# Set column names of the design matrix to match the treatments
colnames(design) <- c("Early_biofilm", "Thin_biofilm", "Mature_biofilm")

# Check the design matrix
design

# Estimate the dispersion for the data
dge <- estimateDisp(dge, design)  # General estimation

# Estimate the robust dispersion to be less sensitive to outliers
dge <- estimateDisp(dge, design, robust=TRUE)

# Get the common, trended, and tagwise dispersions
print(dge$common.dispersion) # 0.03293905
print(dge$trended.dispersion) #0.06886860 0.06780001 0.06565877 0.06560858 more elements ...
print(dge$tagwise.dispersion) # 0.11229441 0.10051499 0.14325575 0.07012374 more elements ...

# Plot the Biological Coefficient of Variation (BCV) versus Average Log CPM
plotBCV(dge) 

#----------------------------------#
# Section 5: Differential expression analysis between all pairs of treatments-----
#----------------------------------#

# Fit the GLM model to the data using the design matrix
fit <- glmQLFit(dge, design)

# Define the contrast matrix
contrast_matrix <- makeContrasts(
  Early_vs_Thin = Early_biofilm - Thin_biofilm,
  Early_vs_Mature = Early_biofilm - Mature_biofilm,
  Thin_vs_Mature = Thin_biofilm - Mature_biofilm,
  levels = design
)

# Perform the differential expression test for all contrasts
res <- glmQLFTest(fit, contrast = contrast_matrix)

# Return the results table
head(res$table)

dim(fit$coefficients) # 5078    3

# Add the gene_symbol to the results
res$table$gene_symbol <- dge$genes$gene_symbol[match(rownames(res$table), dge$genes$gene_id)]

# Replace NA values in the gene_symbol column with "unknown"
res$table$gene_symbol[is.na(res$table$gene_symbol)] <- "unknown"

# Check the first few rows of the result
head(res$table)

## Filter genes with p-value <= 0.05----
de_genes_pvalue_5 <- res$table[res$table$PValue <= 0.05, ]

# Get the total number of differentially expressed genes with p-value <= 0.05
total_de_genes_pvalue_5 <- nrow(de_genes_pvalue_5)

# Print the result
cat("Total number of differentially expressed genes with p-value <= 0.05:", total_de_genes_pvalue_5, "\n") # 3757 

## Apply FDR correction -----
res$table$FDR <- p.adjust(res$table$PValue, method = "BH")

# View the results
topTags(res)

# Extract results
results <- topTags(res, n = Inf)

head(results)

#The total number of differentially expressed genes at 5% FDR 
# Filter the genes based on FDR <= 0.05 (5% significance)
de_genes_fdr_5 <- results$table[results$table$FDR <= 0.05, ]
head(de_genes_fdr_5)

# Get the total number of differentially expressed genes at 5% FDR
total_de_genes_fdr_5 <- nrow(de_genes_fdr_5)

# Print the result
cat("Total number of differentially expressed genes at 5% FDR:", total_de_genes_fdr_5, "\n") #3589 

# Check the distribution of FDR values again
summary(results$table$FDR)

## Extract the list----
# The 'deg_list' contains the filtered DEGs
deg_list_final <- results$table[, c("gene_id", "gene_symbol", "logFC.Early_vs_Thin", "logFC.Early_vs_Mature", "logFC.Thin_vs_Mature", "logCPM", "PValue", "FDR", "F")]

# Show the list of DEGs with F-statistic
head(deg_list_final)

#remove row labels
rownames(deg_list_final) <- NULL

# Save the filtered DEGs to CSV, including the F-statistic
write.csv(deg_list_final, "DEGs_between_treatments.csv", row.names = FALSE)


## Plot heatmap of DEG among the three stages-----
# Calculate logCPM values
logCPM_matrix <- cpm(dge, log = TRUE)

# 1. Sort deg_list_final by absolute logFC for a chosen contrast ("Thin_vs_Mature")
deg_list_final_sorted <- deg_list_final[order(abs(deg_list_final$logFC.Early_vs_Thin), decreasing = TRUE), ]

# 2. Select the top 30 genes based on this metric
top30_deg <- deg_list_final_sorted[1:30, ]
top30_gene_ids <- top30_deg$gene_id

# 3. Check which of these top genes are present in the expression matrix
common_genes <- intersect(top30_gene_ids, rownames(logCPM_matrix))
cat("Number of top genes found in expression data:", length(common_genes), "\n")

# 4. Subset the logCPM matrix to include only the common top genes
heatmap_data <- logCPM_matrix[common_genes, , drop = FALSE]

# Order the heatmap_data to reflect the same order as in top30_deg:
ordered_common <- top30_gene_ids[top30_gene_ids %in% common_genes]
heatmap_data <- heatmap_data[ordered_common, ]

# 5. Prepare annotation for treatments
annotation_col <- data.frame(
  Treatment = sample_metadata$treatments,
  row.names = colnames(heatmap_data)  # Match sample names in heatmap_data
)

# 6. Plot the heatmap using pheatmap with treatment annotations
pheatmap(heatmap_data, 
         scale = "row",  # Scale expression by gene (row)
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         annotation_col = annotation_col,  # Add treatment annotations
         main = "Heatmap of Top 30 Highly Regulated Genes")
#-------------------------------------------------------------#
# Section 6: Genes differentially expressed between the early treatment and both later treatments (Thin and Mature)------
#--------------------------------------------------------------#
# Check the design matrix
print(design)

# Define the contrast: Compare Early biofilm versus the average of Thin and Mature biofilm
contrast_matrix_2 <- makeContrasts(
  Early_vs_Later = Early_biofilm - (Thin_biofilm + Mature_biofilm)/2,
  levels = design
)

print("Contrast matrix_2:")
print(contrast_matrix_2)

# Fit the GLM model using the design matrix
fit2 <- glmQLFit(dge, design)

# Perform the differential expression test for the defined contrast
res2 <- glmQLFTest(fit2, contrast = contrast_matrix_2)

# Return the results table
head(res2$table)

dim(fit2$coefficients) # 5078    3

# Add the gene_symbol to the results
res2$table$gene_symbol <- dge$genes$gene_symbol[match(rownames(res2$table), dge$genes$gene_id)]

# Replace NA values in the gene_symbol column with "unknown"
res2$table$gene_symbol[is.na(res2$table$gene_symbol)] <- "unknown"

# Check the first few rows of the result
head(res2$table)

## Filter genes with p-value <= 0.05----
de_genes2_pvalue_5 <- res2$table[res2$table$PValue <= 0.05, ]

# Get the total number of differentially expressed genes with p-value <= 0.05
total_de_genes2_pvalue_5 <- nrow(de_genes2_pvalue_5)

# Print the result
cat("Total number of differentially expressed genes with p-value <= 0.05:", total_de_genes2_pvalue_5, "\n") # 2990 

## Apply FDR correction -----
res2$table$FDR <- p.adjust(res2$table$PValue, method = "BH")

# Check the distribution of FDR values again
summary(res2$table$FDR)

# View the results
topTags(res2)

# Extract results
results2 <- topTags(res2, n = Inf)

head(results2)

## Extract the list----
# The 'deg_list' contains the filtered DEGs
deg_list_2 <- results2$table[, c("gene_id", "gene_symbol", "logFC", "logCPM", "PValue", "FDR", "F")]

# Show the list of DEGs with F-statistic
head(deg_list_2)

#remove row labels
rownames(deg_list_2) <- NULL

# Save the filtered DEGs to CSV, including the F-statistic
write.csv(deg_list_2, "DEGs_Early_vs_LAter.csv", row.names = FALSE)


## Count number of up- and downregulated genes-----

#The total number of differentially expressed genes at 5% FDR 
# Define thresholds
logFC_threshold <- 1        #  absolute logFC > 1
significance_threshold <- 0.05  # Using FDR 

# Filter genes with FDR <= 0.05 and |logFC| > 1
de_genes2_fdr_5 <- results2$table[
  results2$table$FDR <= 0.05 & abs(results2$table$logFC) > logFC_threshold,
]
head(de_genes2_fdr_5)

# Get the total number of differentially expressed genes at 5% FDR
total_de_genes2_fdr_5 <- nrow(de_genes2_fdr_5)

# Count upregulated genes (logFC > threshold and significant)
upregulated_genes <- results2$table[results2$table$logFC > logFC_threshold & results2$table$FDR < significance_threshold, ]
num_upregulated <- nrow(upregulated_genes)

# Count downregulated genes (logFC < -threshold and significant)
downregulated_genes <- results2$table[results2$table$logFC < -logFC_threshold & results2$table$FDR < significance_threshold, ]
num_downregulated <- nrow(downregulated_genes)

# Print the result
cat("Total number of differentially expressed genes at 5% FDR:", total_de_genes2_fdr_5, "\n") #1021 

# Print the results
cat("Number of upregulated genes:", num_upregulated, "\n") #486 
cat("Number of downregulated genes:", num_downregulated, "\n") #535 


## Create a volcano plot----

# Prepare data for the volcano plot
volcano_data <- results2$table
volcano_data$color <- ifelse(
  volcano_data$logFC > logFC_threshold & volcano_data$FDR < significance_threshold, "Upregulated",
  ifelse(volcano_data$logFC < -logFC_threshold & volcano_data$FDR < significance_threshold, "Downregulated", "Not Significant")
)

#Remove the duplicate "gene_symbol" column
volcano_data <- volcano_data[, !duplicated(colnames(volcano_data))]

# Check the column names again
colnames(volcano_data)

# Plot the volcano plot
ggplot(volcano_data, aes(x = logFC, y = -log10(FDR), color = color)) +
  geom_point(aes(color = color), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: Early vs Later Biofilm",
    x = "Log2 Fold Change (Early vs Later)",
    y = "-Log10 FDR",
    color = "Expression"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the volcano plot
ggsave("volcano_plot_early_vs_later.png", width = 8, height = 6, dpi = 300)


## Plot heatmap DEGs Erly vs Later----

# Calculate logCPM values
logCPM_matrix2 <- cpm(dge, log = TRUE)

# 1. Sort deg_list_2 by absolute logFC 
deg_list_2_sorted <- deg_list_2[order(abs(deg_list_2$logFC), decreasing = TRUE), ]

# 2. Select the top 30 genes based on this metric
top30_deg2 <- deg_list_2_sorted[1:30, ]
top30_gene_ids2 <- top30_deg2$gene_id

# 3. Check which of these top genes are present in the expression matrix
common_genes2 <- intersect(top30_gene_ids2, rownames(logCPM_matrix2))
cat("Number of top genes found in expression data:", length(common_genes2), "\n")

# 4. Subset the logCPM matrix to include only the common top genes
heatmap_data2 <- logCPM_matrix2[common_genes2, , drop = FALSE]

# Order the heatmap_data to reflect the same order as in top30_deg:
ordered_common2 <- top30_gene_ids2[top30_gene_ids2 %in% common_genes2]
heatmap_data2 <- heatmap_data2[ordered_common2, ]

# 5. Prepare annotation for treatments
# Assuming sample_metadata is a dataframe with sample names and treatments
annotation_col <- data.frame(
  Treatment = sample_metadata$treatments,
  row.names = colnames(heatmap_data2)  # Match sample names in heatmap_data2
)

# 6. Plot the heatmap using pheatmap with treatment annotations
pheatmap(heatmap_data2, 
         scale = "row",  
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         annotation_col = annotation_col,  # Add treatment annotations
         main = "Heatmap of Top 30 Highly Regulated Genes (Early vs Later Biofilm)")
