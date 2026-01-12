# Identifying Differentially Expressed Genes Among Three Flor Yeast Growth Stages

In the current study, an experiment mimicking the industrial wine aging, Flor yeast strain was grown, and velum samples were collected at three different growth formations (i.e., early, thin and mature biofilm) to determine the differentially expressed genes at the various formation stages. The goal of this project is to identify the differentially expressed genes (DEGs) among the three velum formation stages and compare the differentially expressed genes at the early (early biofilm) and later stages (thin and mature biofilm). To achieve this, the quality of the RNAseq dataset was determined using fastp, which tests the quality of the sequences and trims the sequences in one step, compared to other tools such as fastqc that require two steps using two different tools. For feature count selection, the multioverlap argument was used, increasing the successfully assigned alignments from an average 15% to 95%. edgeR used for identifying DEGs is a robust method for identifying DEGs, particularly suited for experiments with small sample sizes, as it models biological variation and provides accurate estimates of dispersion (McCarthy et al., 2012). 

## Methods
* Data Acquisition: Three-time points with three biological replicates were used for RNAseq analysis of yeast strain I-329. The short reads (50-bp single-end-read) obtained from Illumina HiSeq 2500 (SRR9027083) and the genome of S. cerevisiae strain I-329 (GenBank PTER00000000) as a reference genome (RG) were used for downstream analysis.

* Quality Control and Trimming: Firstly, fastp was used to perform to test quality control (QC) and trim low-quality bases from the front and tail with a minimum length of 36 bp.

* Alignment and feature count: The trimmed reads were aligned to the pre-indexed reference genome using the STAR aligner. Gene expression quantification was performed using the Rsubread package in R. The featureCounts function, with the allow multi-overlap parameter, was used to count reads mapped to genomic features defined in the GTF annotation file.

* Differential gene expression: Statistical analysis of differential expression was performed using the edgeR package and trimmed mean of M-values (TMM) was used to normalize the transcript counts (Luce et al., 2022). DEGs were identified based on a p-value threshold of ≤ 0.05 and a false discovery rate (FDR) threshold of ≤ 0.05, adjusted using the Benjamini–Hochberg method (Benjamini and Hochberg, 1995). The genes were considered significantly regulated when log2 fold change > 1 and FDR < 0.05. 

<img width="518" height="447" alt="image" src="https://github.com/user-attachments/assets/c1ac3990-a410-41e6-828b-78b4353c4aad" />

Figure 1: Multidimensional Scaling (MDS) plot showing the separation of samples by treatment group (Mature, Thin, Early biofilm).  Samples are colored by treatment groups: Mature biofilm (blue), Thin biofilm (red), and Early biofilm (green) and each point represents a replicate. The x-axis and y-axis represent the leading log-fold changes (logFC) along the first two dimensions, accounting for 59% and 26% of the variation, respectively. 


<img width="477" height="437" alt="image" src="https://github.com/user-attachments/assets/acf63385-7cdf-4f2d-93a9-7efe92861fd7" />

Figure 2: Biological Coefficient of Variation (BCV) plot showing the dispersion estimates for genes across different expression levels. The plot includes tagwise, common, and trended dispersions, illustrating the variability in gene expression. 


<img width="875" height="743" alt="image" src="https://github.com/user-attachments/assets/47c01494-dabc-4457-9357-00a990f9f44f" />

Figure 3: Heatmap of the top 30 highly regulated genes across different treatment groups: Early biofilm, Mature biofilm, and Thin biofilm. The rows represent genes, and the columns represent samples. The color intensity indicates the expression levels (logCPM), with orange representing higher expression and blue representing lower expression.


<img width="448" height="412" alt="image" src="https://github.com/user-attachments/assets/275b154a-3b5b-49fe-bc4b-c6fb1177ee3d" />

Figure 4: Volcano plot comparing gene expression between Early and Later (average of Thin and Mature biofilm) yeast biofilm stages. Points are colored based on significance and direction of change: Upregulated genes (red) have logFC > 1 and FDR < 0.05, Downregulated genes (blue) have logFC < -1 and FDR < 0.05, and Not Significant genes (gray) do not meet these thresholds. 

<img width="840" height="713" alt="image" src="https://github.com/user-attachments/assets/938db803-c441-4339-b372-55985f86d737" />

Figure 5: Heatmap of the top 30 highly regulated genes in the comparison between Early and Later (average of Thin and Mature biofilm) yeast biofilm stages. The rows represent genes, and the columns represent samples grouped by treatment. The color intensity indicates the expression levels (logCPM), with orange representing higher expression and blue representing lower expression. 







