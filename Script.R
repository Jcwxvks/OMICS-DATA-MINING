library(data.table)
library(tximport)
library(edgeR)
library(reshape2)
library(jsonlite)
library(ggrepel)
library(edgeR)
library(ggplot2)
library(writexl)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(topGO)
library(org.Hs.eg.db)

# read the phenotype data
pheno <- fread("phenotype.txt")
setkey(pheno, "sample")

# get the annotations using one of the runs
annotations <- read.csv(sprintf("kallisto_run_%s/abundance.tsv", pheno$sample[[1]]), sep="\t", header=T, quote="", as.is=T)[1]
annotations$TXNAME <- sapply(strsplit(annotations[[1]], "|", fixed=T), function(x) x[1])
annotations$GENEID <- sapply(strsplit(annotations[[1]], "|", fixed=T), function(x) x[2])
annotations$gene_name <- sapply(strsplit(annotations[[1]], "|", fixed=T), function(x) x[6])
annotations$gene_type <- sapply(strsplit(annotations[[1]], "|", fixed=T), function(x) x[8])
annotations <- annotations[annotations$gene_type=="protein_coding", ]

# import the transcript level data
pnames <- sprintf("kallisto_run_%s/abundance.tsv", pheno$sample)
names(pnames) <- pheno$sample
txi <- tximport(pnames, type="kallisto", tx2gene=annotations[, c("TXNAME", "GENEID")], ignoreAfterBar=T)

# get the mapping rate
pnames <- sprintf("kallisto_run_%s/run_info.json", pheno$sample)
mapping_data <- cbind(sample=gsub("^kallisto_run_", "", pnames), do.call(rbind, lapply(lapply(pnames, fromJSON), as.data.frame)))
mapping_data$call <- NULL

# prepare the data for creating an object for use with edgeR
cts <- txi$counts
normMat <- txi$length

# obtaining per-observation scaling factors for length, adjusted to avoid changing the magnitude of the counts
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# computing effective library sizes from scaled counts, to account for composition biases between samples
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# combining effective library sizes with the length factors, and calculating offsets for a log-link GLM
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# creating a DGEList object for use in edgeR.
e <- DGEList(cts)
e <- scaleOffset(e, normMat)
e$samples$cell_type <- factor(pheno[rownames(e$samples), ]$cell_type, levels=c("Basophils", "DN_GMPs"))
e$samples$sample <- rownames(e$samples)
fd <- unique(annotations[, c("GENEID", "gene_name", "gene_type")])
colnames(fd) <- c("gene_id", "gene_name", "gene_type")
rownames(fd) <- fd$gene_id
fd$length <- apply(txi$length, 1, mean)
e$genes <- fd[rownames(e$counts), ]

# save the object
save(e, file="rna_seq_data.RData")

# generate the log RPKM
e_rpkm <- rpkm(e, log=T)
e_rpkm <- melt(e_rpkm)
colnames(e_rpkm) <- c("gene_id", "sample", "rpkm")
e_rpkm <- merge(e_rpkm, e$samples[, c("sample", "cell_type")], by="sample")
e_rpkm <- merge(e_rpkm, e$genes, by="gene_id")
counts_data <- melt(txi$counts)
colnames(counts_data) <- c("gene_id", "sample", "count")
e_rpkm <- merge(e_rpkm, counts_data, by=c("gene_id", "sample"))

# boxplot
pdf("boxplots.pdf", width=8, height=6)
p <- ggplot(e_rpkm, aes(x=sample, y=count)) + geom_boxplot() + theme_bw() + ggtitle("Box plot of read counts") + ylab("read counts") + scale_y_log10()
print(p)
p <- ggplot(e_rpkm, aes(x=sample, y=rpkm)) + geom_boxplot() + theme_bw() + ggtitle("Box plot of log2 RPKM") + ylab("log2 RPKM")
print(p)
dev.off()

# run the PCA
d1 <- acast(e_rpkm, sample ~ gene_id, value.var="rpkm")
fit <- prcomp(d1)
pca_pcs <- cbind(e$samples, as.data.frame(fit$x))
pca_var <- data.frame(pc=1:length(fit$sdev), variance=(fit$sdev^2)/sum(fit$sdev^2))
pca_loadings <- cbind(parameter=rownames(fit$rotation), as.data.frame(fit$rotation))

# plot the PCA
pdf("pca.pdf", width=8, height=6)
p <- ggplot(pca_var, aes(x=pc, y=variance*100)) + geom_col() + theme_bw() + ggtitle("Variance accounted for")
print(p)
p <- ggplot(pca_pcs, aes(x=PC1, y=PC2, col=cell_type, label=sample)) + geom_point() + geom_text_repel(col="black") + theme_bw() + ggtitle("PCA plot") + xlab(sprintf("PC1 (%0.3f%%)", pca_var[1, 2]*100)) + ylab(sprintf("PC2 (%0.3f%%)", pca_var[2, 2]*100))
print(p)
dev.off()

# plot the number of reads per sample
pdf("reads_bar_chart.pdf", width=12, height=6)
p <- ggplot(mapping_data, aes(x=sample, y=n_processed)) + geom_col() + theme_bw() + ggtitle("Number of reads per sample")
print(p)
dev.off()

# plot the mapping rate per sample
pdf("mapping_rates_bar_chart.pdf", width=12, height=6)
p <- ggplot(mapping_data, aes(x=sample, y=p_pseudoaligned)) + geom_col() + theme_bw() + ggtitle("Mapping rate per sample")
print(p)
dev.off()

# save the data
write_xlsx(list(
  Data=e_rpkm,
  PC=pca_pcs,
  Variance=pca_var,
  Loadings=pca_loadings
), "qc_data.xlsx")

# run edgeR DEG analysis
design <- model.matrix(~e$samples$cell_type)
e <- estimateDisp(e, design=design)
fit <- glmQLFit(e, design)
results <- as.data.frame(topTags(glmQLFTest(fit), n=Inf, adjust.method="BH"))
results$gene_id <- mapIds(org.Hs.eg.db, keys=results$gene_name, keytype="SYMBOL", column="ENTREZID", multiVals="first")
results$rank <- rank(results$PValue)
results <- results[order(results$PValue), ]

# save the results to an Excel file
write_xlsx(list(DEGs=results), "edger_results.xlsx")

# volcano plot
results$state <- ifelse(results$FDR<0.05, "sig", "not_sig")
results$label <- ifelse(results$rank <= 25, results$gene_name, NA_character_)
pdf("edger_volcano_plot.pdf", width=12, height=10)
p <- ggplot(results, aes(x=logFC, y=-log10(PValue), col=-log10(FDR), label=label)) + geom_point() + scale_color_gradient2(low="grey", mid="grey", high="red", midpoint=-log10(0.05)) + geom_vline(xintercept=0, size=0.5, linetype="dotted") + geom_hline(yintercept=0, size=0.5, linetype="dotted") + theme_bw() + geom_text_repel(col="black") + ggtitle("Volcano plot")
print(p)
dev.off()

# heat map
d2 <- acast(e_rpkm[e_rpkm$gene_id %in% results[results$rank<=10, ]$gene_id, ], gene_name ~ sample, value.var="rpkm")
d2 <- t(scale(t(d2)))
pdf("edger_heatmap.pdf", height=4 + (nrow(d2)*0.1), width=3 + (ncol(d2)*0.2))
p <- Heatmap(d2, row_names_gp=gpar(fontsize=6), column_names_gp=gpar(fontsize=6), clustering_distance_columns="euclidean", clustering_distance_rows="euclidean", column_dend_side="top", column_names_side="top", column_dend_height=unit(30, "mm"), row_dend_width=unit(30, "mm"), row_names_side="left", col=colorRamp2(c(min(d2), 0, max(d2)), c("blue", "white", "red")))
print(p)
dev.off()

# topGO gene ontology enrichment
all_genes <- ifelse(results$FDR<0.05, 1, 0)
names(all_genes) <- results$entrez_id
all_genes <- all_genes[!is.na(results$entrez_id)]
d1 <- new("topGOdata", ontology="BP", allGenes=all_genes, geneSel=function(x) { return(x==1) }, annot=annFUN.org, nodeSize=10, mapping="org.Hs.eg.db", ID="entrez")
weight01 <- runTest(d1, algorithm = "weight01", statistic = "fisher")
go_results <- GenTable(d1, weight01=weight01, orderBy="weight01", topNodes=length(score(weight01)), numChar=10000)
go_results$weight01 <- as.numeric(gsub("<\\s+", "", go_results$weight01))
go_results <- rename(go_results, c(GO.ID="goid", Term="description", Annotated="annotated", Significant="significant", Expected="expected", weight01="pvalue"))
go_results$fdr <- p.adjust(go_results$pvalue, "BH")

# get the gene mapping data
go_mapping <- data.table()
if (nrow(go_results[go_results$pvalue<0.05, ]) > 0) {
  go_mapping <- do.call(rbind, lapply(go_results[go_results$pvalue<0.05, ]$goid, function(x) {
    gnames <- genesInTerm(d1, x)[[1]]
    return(data.frame(goid=x, gene_id=gnames[gnames %in% names(all_genes[all_genes==1])]))
  }))
  go_mapping <- merge(go_results[, c("goid", "description")], go_mapping, by="goid")
  go_mapping$gene_name <- mapIds(org.Hs.eg.db, keys=go_mapping$gene_id, keytype="geneID", column="SYMBOL")
}

# save the results to an Excel file
write_xlsx(list(
  GO=go_results,
  Mapping=go_mapping
), "topgo_results.xlsx")