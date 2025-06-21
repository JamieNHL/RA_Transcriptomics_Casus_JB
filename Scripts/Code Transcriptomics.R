setwd("C:/Users/jamie/Downloads/Transcriptomics Project")
getwd()
library(Rsubread)

#Indexing
buildindex(
  basename = 'ref_human',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE
)

#Aligning
align.con1 <- align(index = "ref_human", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2 = "SRR4785819_2_subset40k.fastq", output_file = "con1.BAM")
list.files(pattern = "SRR4785819")
align.con2 <- align(index = "ref_human", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2 = "SRR4785820_2_subset40k.fastq", output_file = "con2.BAM")
align.con3 <- align(index = "ref_human", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2 = "SRR4785828_2_subset40k.fastq", output_file = "con3.BAM")
align.con4 <- align(index = "ref_human", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2 = "SRR4785831_2_subset40k.fastq", output_file = "con4.BAM")
align.ra1  <- align(index = "ref_human", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2 = "SRR4785979_2_subset40k.fastq", output_file = "ra1.BAM")
align.ra2  <- align(index = "ref_human", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2 = "SRR4785980_2_subset40k.fastq", output_file = "ra2.BAM")
align.ra3  <- align(index = "ref_human", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2 = "SRR4785986_2_subset40k.fastq", output_file = "ra3.BAM")
align.ra4  <- align(index = "ref_human", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2 = "SRR4785988_2_subset40k.fastq", output_file = "ra4.BAM")

#Sorting
library(Rsamtools)
samples <- c('con1', 'con2', 'con3', 'con4', 'ra1', 'ra2', 'ra3', 'ra4')
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})

lapply(samples, function(s) {indexBam(file = paste0(s, '.sorted.bam'))
})

#Count Matrix
library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)
allsamples <- c(
  "con1.sorted.bam", "con2.sorted.bam", "con3.sorted.bam", "con4.sorted.bam",
  "ra1.sorted.bam", "ra2.sorted.bam", "ra3.sorted.bam", "ra4.sorted.bam"
)
count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "GRCh38.p14_genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)
str(count_matrix)
head(count_matrix$counts)
counts <- read.table("count_matrix.txt", 
                     header = TRUE, 
                     row.names = 1, 
                     sep = "\t")
head(counts)
colnames(counts)
treatment <- c("control", "control", "control", "control",
               "RA", "RA", "RA", "RA")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- colnames(counts)

#Data Analysis
BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
library(DESeq2)
library(KEGGREST)

#DESeq Data Analysis
dds <- DESeqDataSetFromMatrix(round(counts),
                              colData = treatment_table,
                              design = ~ treatment)
summary(as.vector(as.matrix(counts)))
dds <- DESeq(dds)
resultaten <- results(dds)
write.table(resultaten, file = 'DESeq2_Resultaten.csv', row.names = TRUE, col.names = TRUE)

#Differentiële expressie
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]
head(laagste_p_waarde)

#Volcano Plot
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)
if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}
library(pathview)

#KEGG Analysis
res_df <- as.data.frame(resultaten) <-- #DESeq2-Resultaat
significant_genen <- res_df[which(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1), ]
head(significant_genen)
nrow(significant_genen)
str(significant_genen)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(clusterProfiler)
library(org.Hs.eg.db)
head(rownames(significant_genen))
entrez_ids <- bitr(
  rownames(significant_genen),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
kegg_enrich <- enrichKEGG(
  gene         = entrez_ids$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.05
)
head(kegg_enrich)
kegg_df <- as.data.frame(kegg_enrich)
write.csv(kegg_df, file = "KEGG_enrichment.csv", row.names = FALSE)
dotplot(kegg_enrich, showCategory = 10)
ggsave("kegg_dotplot.png", plot = last_plot(), width = 10, height = 8, dpi = 300)

#GO Analysis
if (!requireNamespace("goseq", quietly = TRUE)) {
  BiocManager::install("goseq")
}
if (!requireNamespace("GO.db", quietly = TRUE)) {
  BiocManager::install("GO.db")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(goseq)
library(GO.db)
library(tidyverse)
load("Robjects/Annotated_Results_LvV.RData")
is_significant <- as.integer(resultaten$padj < 0.05 & abs(resultaten$log2FoldChange) > 1)
names(is_significant) <- rownames(resultaten)
supportedOrganisms() %>% filter(str_detect(Genome, "hg"))
pwf <- nullp(
  DEgenes = is_significant,
  genome = "hg19",
  id = "geneSymbol"
)
goResults <- goseq(pwf, "hg19", "geneSymbol", test.cats = c("GO:BP"))
write.csv(goResults, file = "GO_results.csv", row.names = FALSE)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)
library(AnnotationDbi)
gene_symbols <- rownames(resultaten)
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
is_significant_entrez <- as.integer(resultaten$padj < 0.05 & abs(resultaten$log2FoldChange) > 1)
names(is_significant_entrez) <- entrez_ids
is_significant_entrez <- is_significant_entrez[!is.na(names(is_significant_entrez))]
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
}
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gene_lengths <- transcriptLengths(txdb) %>%
  as_tibble() %>%
  group_by(gene_id) %>%
  summarise(length = sum(tx_len)) %>%
  distinct(gene_id, .keep_all = TRUE)

length_vector <- gene_lengths$length
names(length_vector) <- gene_lengths$gene_id

length_vector <- length_vector[names(is_significant_entrez)]
valid_genes <- !is.na(length_vector)
length_vector <- length_vector[valid_genes]
is_significant_entrez <- is_significant_entrez[valid_genes]

pwf <- nullp(DEgenes = is_significant_entrez, bias.data = length_vector)

sum(is.na(is_significant_entrez))
sum(is.na(length_vector))

identical(names(is_significant_entrez), names(length_vector))

length(is_significant_entrez)
length(length_vector)

valid <- !is.na(is_significant_entrez)

is_significant_entrez <- is_significant_entrez[valid]
length_vector <- length_vector[valid]

sum(is.na(is_significant_entrez))
sum(is.na(length_vector))

pwf <- nullp(DEgenes = is_significant_entrez, bias.data = length_vector)

goResults <- goseq(pwf, "hg19", "entrez", test.cats = c("GO:BP"))

library(org.Hs.eg.db)

GO_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = names(is_significant_entrez),
                                        columns = c("GO"),
                                        keytype = "ENTREZID")

GO_annotations <- GO_annotations[GO_annotations$ONTOLOGY == "BP", ]

goResults <- goseq(pwf,
                   gene2cat = split(GO_annotations$GO, GO_annotations$ENTREZID))
library(ggplot2)

goResults %>%
  top_n(10, wt = -over_represented_pvalue) %>%
  mutate(hitsPerc = numDEInCat * 100 / numInCat) %>%
  ggplot(aes(x = hitsPerc,
             y = reorder(term, hitsPerc),
             colour = over_represented_pvalue,
             size = numDEInCat)) +
  geom_point() +
  expand_limits(x = 0) +
  labs(x = "Hits (%)",
       y = "GO term",
       colour = "p-value",
       size = "Differentially expressed")
ggsave("GO_plot.png", width = 10, height = 6, dpi = 300)
