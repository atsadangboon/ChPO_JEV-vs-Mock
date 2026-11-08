library("Rsubread")
library("limma")
library("edgeR")

directory <-
  "/yourdirectory"
setwd(directory)

gtf = "hg38.ncbiRefSeq.gtf"
bams = list.files()
bams = grep("*.bam$", bams, value=TRUE)
Fc <- featureCounts(
  bams,
  annot.ext           = gtf,
  isGTFAnnotationFile = TRUE,
  isPairedEnd         = TRUE,
  GTF.featureType     = "exon",
  GTF.attrType        = "gene_id",
  useMetaFeatures     = FALSE,
  allowMultiOverlap   = TRUE,
  nthreads            = 4,
  
)

#write table for collect Fc output
write.table(x=data.frame(Fc$annotation[,c("GeneID","Length")], Fc$counts, stringsAsFactors=FALSE), file="rawcount.csv", quote=FALSE, sep="\t", row.names=FALSE,)


df <- read.table("rawcount.csv",header = TRUE, sep = "\t")
head(df) #check table file
colnames(df) [3] <- "Counts" #rename colname 

#rename colname 
colnames(df) [3] <- "J1"
colnames(df) [4] <- "J2"
colnames(df) [5] <- "J3"
colnames(df) [6] <- "M1"
colnames(df) [7] <- "M2"
colnames(df) [8] <- "M3"

# Exclude columns 9 to 11
df <- df[, -c(9:11)]

# Check the updated data frame
head(df) #check output

#aggregate duplicate data
df1 <- aggregate(df[,3:8],by = list(df$GeneID), FUN=sum)
head(df1,2) #check output

#set GeneID as row.name
sum_counts <- df1[,-1]
rownames(sum_counts) <- df1[,1]
head(sum_counts,2) #check output

keep <- rowSums(sum_counts>10) >=2
sum_counts <- sum_counts[keep,]

library(DESeq2) #to run package, make sure you recall it
count_matrix <- sum_counts #count matrix from read counts
head(count_matrix) #To check count_matrix

#To define which samples are in each column
sample = c("J1", "J2", "J3", "M1", "M2", "M3")
condition = c("J", "J", "J", "M", "M", "M")
coldata <- data.frame(sample, condition, row.names = "sample")
coldata$condition <- as.factor(coldata$condition)
condition
head(coldata)

#Check matching rownames and colnames (should be true)
all(rownames(coldata) %in% colnames(count_matrix)) #check TRUE
all(rownames(coldata) == colnames(count_matrix)) #check TRUE

#Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                              colData = coldata, 
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "M") #set M as control

#Run DESeq2
dds<-DESeq(dds)

#Get results
resultsNames(dds) #see all comparisons
res <- results(dds)
head(res) #Inspect results

summary(results(dds, alpha = 0.05)) #summary of DEG

# Filter results for padj < 0.05
significant_res <- res[which(res$padj < 0.05), ]
head(significant_res)

# Count the number of significant genes
num_significant_genes <- nrow(significant_res)
cat("Number of significant genes (padj < 0.05):", num_significant_genes, "\n")

# Filter genes with |log2FC| > 1 and padj < 0.05
filtered_res <- res[which(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1.0), ]
head(filtered_res)
num_filtered_genes <- nrow(filtered_res)
cat("Number of genes with |log2FC| > 1 and padj < 0.05:", num_filtered_genes, "\n")


normalized_counts <- counts(dds, normalized=TRUE) #Normalized counts
head(normalized_counts,100) #show result


library(dplyr)

volcano_data <- as.data.frame(res)
volcano_data$gene <- rownames(volcano_data)
volcano_data <- volcano_data %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
    padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

highlight_genes <- c("CXCL6", "CXCL10", "IFNB1")

ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
  geom_point(alpha=0.6, size=3) +
  scale_color_manual(values=c("Upregulated"="#48C9B0", "Downregulated"="#E74C3C", "Not Significant"="gray")) +
  geom_text(data=filter(volcano_data, gene %in% highlight_genes),
            aes(label=gene), nudge_y=3, color="black", size=3) +
  geom_segment(data=filter(volcano_data, gene %in% highlight_genes),
               aes(x=log2FoldChange + 0.2, xend=log2FoldChange,
                   y=-log10(padj) + 3, yend=-log10(padj)),
               color="black", linewidth=0.5, lineend="round") +
  coord_cartesian(xlim=c(-10, 10), ylim=c(0, 45)) +
  theme_minimal(base_size = 14) +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 Adjusted p-value") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank())


library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggforce)

# Heatmap of DEGs (padj < 0.05 & |log2FC| > 1)
de_genes <- rownames(res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ])
vsd <- vst(dds, blind = FALSE)
# Now you can use assay() on vsd
assay(vsd)
vsd_mat <- assay(vsd)[de_genes, ]


# z-score scaling by row
vsd_mat_scaled <- t(scale(t(vsd_mat)))
pheatmap(vsd_mat_scaled,
         annotation_col = coldata,
         show_rownames = FALSE,
         fontsize_col = 12,
         main = "Heatmap of DEGs \n(padj < 0.05 & |log2FC| > 1)")

library(ggplot2)
library(ggrepel)

# Define genes to highlight
highlight_genes <- c("ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX", "BCAP31", 
                     "BCL10", "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN", "BID", "BIK", "BIRC3", 
                     "BMF", "BMP2", "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP1", "CASP2", "CASP3", 
                     "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CAV1", "CCNA1", "CCND1", "CCND2", 
                     "CD14", "CD2", "CD38", "CD44", "CD69", "CDC25B", "CDK2", "CDKN1A", "CDKN1B", 
                     "CFLAR", "CLU", "CREBBP", "CTH", "CTNNB1", "CYLD", "DAP", "DAP3", "DCN", "DDIT3", 
                     "DFFA", "DIABLO", "DNAJA1", "DNAJC3", "DNM1L", "DPYD", "EBP", "EGR3", "EMP1", 
                     "ENO2", "ERBB2", "ERBB3", "EREG", "ETF1", "F2", "F2R", "FAS", "FASLG", "FDXR", 
                     "FEZ1", "GADD45A", "GADD45B", "GCH1", "GNA15", "GPX1", "GPX3", "GPX4", "GSN", 
                     "GSR", "GSTM1", "GUCY2D", "H1-0", "HGF", "HMGB2", "HMOX1", "HSPB1", "IER3", 
                     "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18", "IL1A", "IL1B", "IL6", 
                     "IRF1", "ISG20", "JUN", "KRT18", "LEF1", "LGALS3", "LMNA", "PLPPR4", "LUM", 
                     "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4", "PDGFRB", 
                     "PEA15", "PLAT", "PLCB2", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1", "PSEN1", 
                     "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2", "RNASEL", "ROCK1", 
                     "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1", "SQSTM1", 
                     "TAP1", "TGFB2", "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", 
                     "TNFSF10", "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP")

# Subset and mark highlight genes
highlight_df <- volcano_data[rownames(volcano_data) %in% highlight_genes, ]
volcano_data$highlight <- ifelse(rownames(volcano_data) %in% highlight_genes, "yes", "no")

p <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  # Main volcano points
  geom_point(aes(color = sig), alpha = 0.6, size = 3) +
  
  # Highlight gene circles
  geom_point(
    data = highlight_df,
    shape = 21, fill = NA, color = "black", stroke = 2, size = 3
  ) +
  
  # Axis cutoffs
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  # Color scale
  scale_color_manual(values = c(
    "Upregulated" = "#d73027",
    "Downregulated" = "#4575b4",
    "Not Significant" = "gray"
  )) +
  
  coord_cartesian(xlim = c(-12, 12), ylim = c(0, 100)) +
  labs(
    title = "Apoptosis",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

print(p)


#GSEA analysis
# ─── 1) load all required libraries ──────────────────────────────────────────
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(msigdbr)
library(AnnotationDbi)
library(openxlsx)
library(dplyr)
library(DOSE)
library(tibble)
library(tidyr)

# ─── 2) read your DESeq2 output ──────────────────────────────────────────────
setwd("/Volumes/AtsadangSSD/Pin/MayDE/JM/")
df <- read.table("DEGs_JM_NA.csv", header = TRUE, row.names = 1)

# ─── 3) map SYMBOL → ENTREZID and merge with log2FC ──────────────────────────
gene_ids <- bitr(rownames(df),
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)
df2 <- merge(gene_ids, df,
             by.x = "SYMBOL",
             by.y = "row.names")

# ─── 4) build & sort the named vector (Entrez → log2FC) ──────────────────────
ranked_genes <- sort(
  setNames(df2$log2FoldChange, df2$ENTREZID),
  decreasing = TRUE
)

# ─── 5) pull Hallmark TERM2GENE from msigdbr ─────────────────────────────────
msig <- msigdbr(species    = "Homo sapiens", 
                collection = "H")
msig$ENTREZID <- mapIds(org.Hs.eg.db,
                        keys      = msig$gene_symbol,
                        column    = "ENTREZID",
                        keytype   = "SYMBOL",
                        multiVals = "first")
term2gene <- subset(msig, !is.na(ENTREZID))[, c("gs_name", "ENTREZID")]

# ─── 6) run GSEA ONCE (keep everything) ───────────────────────────────────────
gseaCp <- GSEA(
  geneList      = ranked_genes,
  TERM2GENE     = term2gene,
  pAdjustMethod = "BH",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,      # <-- keep all pathways
  verbose       = FALSE
)

# ─── 7) enrichment plots ─────────────────────────────────────────────────────

## 7a) top pathway
topID <- gseaCp@result$ID[1]
gseaplot2(gseaCp,
          geneSetID = topID,
          title     = topID)

## 7b) a specific pathway
gseaplot2(gseaCp,
          geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
          title     = "TNFα Signaling via NF-κB",
          color     = "green",
          base_size = 12)

## 7c) loop over all significant (adjusted-p < 0.05) pathways and save PNGs
sigIDs <- gseaCp@result$ID[gseaCp@result$p.adjust < 0.05]
for(id in sigIDs) {
  p <- gseaplot2(gseaCp,
                 geneSetID = id,
                 title     = id,
                 color     = "green",
                 base_size = 12)
  ggsave(filename = paste0(id, ".png"),
         plot     = p,
         width    = 5.5,
         height   = 5)
}

## 7d) dotplot of top 10 by NES (all p‐values)
dotplot(gseaCp, showCategory = 10) +
  ggtitle("GSEA Dotplot (Hallmark)")

## 7e) dotplot of only significant pathways
#    — simply filter the single GSEAResult object
gseaCpSig <- gseaCp
gseaCpSig@result <- subset(gseaCp@result, p.adjust < 0.05)

dotplot(gseaCpSig, showCategory = 11) +
  ggtitle("GSEA J vs M (padj ≤ 0.05)") +
  theme(
    axis.text.y = element_text(size = 10, family = "Helvetica"),
    plot.title  = element_text(size = 12, face = "bold", hjust = 0.5)
  )

#dotplot with NES in x axis
library(dplyr)
library(ggplot2)

# 1) pull your GSEA result into a data.frame
res_df <- as.data.frame(gseaCpSig@result)

# 2) filter, sort, compute ratio, take top 11, reverse factor levels
res_top11 <- res_df %>%
  filter(p.adjust <= 0.05) %>%
  arrange(desc(NES)) %>%
  mutate(
    GeneCount   = lengths(strsplit(core_enrichment, "/")),
    GeneRatio   = GeneCount / setSize
  ) %>%
  slice_head(n = 11) %>%
  mutate(
    Description = factor(Description, levels = rev(Description))
  )

# 3) plot with panel border
# 1) Build and store your plot
p <- ggplot(res_top11, aes(
  x     = NES,
  y     = Description,
  size  = GeneRatio,
  color = p.adjust
)) +
  geom_point() +
  scale_color_gradient(
    low  = "#E13423",
    high = "#095BAE",
    name = "padj"
  ) +
  labs(
    title = "GSEA J vs M (padj ≤ 0.05)",
    x     = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 8, family = "Helvetica"),
    axis.text.x  = element_text(size = 10),
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

print(p)
