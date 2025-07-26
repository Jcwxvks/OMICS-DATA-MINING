# Cell Differentiation Analysis: GMP to Basophil Development

This repository provides an analysis pipeline to investigate the transcriptional changes during the differentiation of granulocyte-monocyte progenitors (GMP) into mature basophils in mice, using the public GEO dataset [GSE132122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132122).

## Project Overview

Cell differentiation is a crucial process by which progenitor cells generate mature cells with specialized functions. This project focuses on understanding the molecular changes involved in the transition from GMP to basophils.

## Repository Contents

- **Kallisto Outputs**  
  Quantification results from running Kallisto on the RNA-Seq samples.

- **Sample Metadata**  
  A tab-separated file containing phenotype information for 6 samples:
  - 3 GMP samples
  - 3 mouse basophil samples

- **`Script.R`**  
  An R script that performs the complete downstream analysis:
  - Reads gene expression data and phenotype information.
  - Creates a `DGEList` object and saves it as `rna_seq.RData`.
  - Generates a boxplot of log2 RPKM values for each sample (`boxplots.pdf`).
  - Performs PCA and plots the first two principal components colored by cell type (`pca_plot.pdf`).
  - Plots the number of reads per sample (`reads_bar_chart.pdf`).
  - Plots the mapping rates per sample (`mapping_rates_bar_chart.pdf`).
  - Runs differential expression analysis using **edgeR** to compare GMP vs. basophils.
  - Saves DEG results to `edger_results.xlsx`.
  - Generates a volcano plot of DEG results (`edger_volcan_plot.pdf`).
  - Creates a heatmap of the top 10 DEGs (`edger_heatmap.pdf`).
  - Performs Gene Ontology enrichment analysis for significant DEGs.
  - Saves GO enrichment results to `go_enrichment_results.xlsx`.

## How to Use

1. **Clone the repository:**
   ```bash
   git clone https://github.com/<your-username>/<your-repo>.git
2. **Install required R packages**, for example,
   ```R
   install.packages("BiocManager")
   BiocManager::install("edgeR")
   install.packages("ggplot2")
   install.packages("readxl")
   # Add other packages as needed
3. **Run the analysis:**
   ```R
   source("Script.R")
All intermediate and final results will be saved in the working directory.

## References

- **GEO Dataset: [GSE132122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132122)**
- **Tools: [Kallisto](https://pachterlab.github.io/kallisto/), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)**

## License

This project is open for academic and educational purposes. Please cite the original GEO dataset if you use this analysis in your work.
