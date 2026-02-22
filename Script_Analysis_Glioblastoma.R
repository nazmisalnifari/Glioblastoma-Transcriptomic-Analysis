# ==========================================
# 1. PANGGIL LIBRARY
# ==========================================
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133plus2.db)
library(clusterProfiler)
library(org.Hs.eg.db)

# ==========================================
# 2. DOWNLOAD & PRE-PROCESSING DATA
# ==========================================
message("Sedang mendownload data GSE50161, mohon tunggu sebentar...")
gset <- getGEO("GSE50161", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# [PENGAMANAN: HANYA AMBIL GLIOBLASTOMA & NORMAL]
group_info <- pData(gset)[["source_name_ch1"]]
idx_gbm <- grep("glioblastoma", group_info, ignore.case = TRUE)
idx_normal <- grep("normal", group_info, ignore.case = TRUE)

# Potong dataset agar bersih dari kanker jenis lain
gset_sub <- gset[ , c(idx_gbm, idx_normal)]
ex <- exprs(gset_sub)

# Log2 Transformasi 
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

# ==========================================
# 2.5 EXPLORATORY DATA ANALYSIS (EDA)
# ==========================================
# Pastikan library umap terinstall untuk visualisasi UMAP
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}
library(umap)

# --- A. BOXPLOT ---
# Untuk melihat apakah distribusi data antar sampel sudah sejajar/normal
group_colors <- ifelse(gset_sub$group == "Glioblastoma", "salmon", "turquoise")
boxplot(ex, col = group_colors, las = 2, outline = FALSE,
        main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
        ylab = "Expression Value (log2)", cex.axis = 0.6)
legend("topright", legend = levels(gset_sub$group), fill = c("salmon", "turquoise"), cex = 0.8)

# --- B. DENSITY PLOT ---
# Untuk melihat kurva sebaran nilai ekspresi secara global
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset_sub$group, each = nrow(ex))
)
plot_density <- ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(title = "Distribusi Nilai Ekspresi Gen (Glioblastoma vs Normal)",
       x = "Expression Value (log2)", y = "Density")
print(plot_density)

# --- C. UMAP PLOT ---
# Untuk melihat apakah sampel Glioblastoma dan Normal mengelompok terpisah
umap_input <- t(ex)
umap_result <- umap(umap_input)
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset_sub$group
)
plot_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP Plot: Glioblastoma vs Normal", x = "UMAP 1", y = "UMAP 2")
print(plot_umap)
# =========================================
# 3. STATISTIK DEG (LIMMA)
# ==========================================
# Beri nama grup yang jelas dan anti-tertukar
gset_sub$group <- ifelse(grepl("normal", pData(gset_sub)[["source_name_ch1"]], ignore.case=TRUE), "Normal", "Glioblastoma")
gset_sub$group <- factor(gset_sub$group, levels = c("Glioblastoma", "Normal"))

design <- model.matrix(~0 + gset_sub$group)
colnames(design) <- levels(gset_sub$group)

# Buat matriks perbandingan murni: Glioblastoma vs Normal
fit <- lmFit(ex, design)
contrast_matrix <- makeContrasts("Glioblastoma - Normal", levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

topTableResults <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf, p.value = 0.01)

# ==========================================
# 4. ANOTASI GEN & VOLCANO PLOT
# ==========================================
probe_ids <- rownames(topTableResults)
gene_annotation <- AnnotationDbi::select(hgu133plus2.db, keys = probe_ids, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")
topTableResults$PROBEID <- rownames(topTableResults)
topTableResults <- merge(topTableResults, gene_annotation, by = "PROBEID", all.x = TRUE)

volcano_data <- data.frame(logFC = topTableResults$logFC, adj.P.Val = topTableResults$adj.P.Val, Gene = topTableResults$SYMBOL)
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01] <- "DOWN"

plot_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) + scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") + geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() + ggtitle("Volcano Plot: Glioblastoma vs Normal")
print(plot_volcano)

# ==========================================
# 5. HEATMAP (Top 50 DEGs)
# ==========================================
topTableResults <- topTableResults[order(topTableResults$adj.P.Val), ]
top50 <- head(topTableResults, 50)
mat_heatmap <- ex[top50$PROBEID, ]
rownames(mat_heatmap) <- ifelse(is.na(top50$SYMBOL) | top50$SYMBOL == "", top50$PROBEID, top50$SYMBOL)
mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]
mat_heatmap <- mat_heatmap[apply(mat_heatmap, 1, var) > 0, ]

annotation_col <- data.frame(Group = gset_sub$group)
rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(mat_heatmap, scale = "row", annotation_col = annotation_col, show_colnames = FALSE, show_rownames = TRUE, fontsize_row = 7, main = "Top 50 DEGs in Glioblastoma")

# ==========================================
# 6. ENRICHMENT ANALYSIS (GO & KEGG)
# ==========================================
significant_genes <- topTableResults$SYMBOL[volcano_data$status %in% c("UP", "DOWN")]
significant_genes <- na.omit(significant_genes)
entrez_ids <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

ego <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
plot_go <- dotplot(ego, title = "Gene Ontology (GO) - Glioblastoma")
print(plot_go)

ekegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
plot_kegg <- dotplot(ekegg, title = "KEGG Pathway - Glioblastoma")
print(plot_kegg)

# Export data ke CSV
write.csv(topTableResults, "GSE50161_Glioblastoma_vs_Normal_Results.csv", row.names = FALSE)
message("Analisis Selesai! Keempat grafik sudah siap di panel 'Plots'.")

# --- PERBAIKAN MUTLAK: Mendefinisikan ulang grup agar R tidak bingung ---
gset_sub$group <- factor(ifelse(grepl("normal", pData(gset_sub)[["source_name_ch1"]], ignore.case=TRUE), "Normal", "Glioblastoma"), levels = c("Glioblastoma", "Normal"))

# --- 1. BOXPLOT WARNA-WARNI ---
group_colors <- ifelse(gset_sub$group == "Glioblastoma", "salmon", "turquoise")
boxplot(ex, col = group_colors, las = 2, outline = FALSE,
        main = "Boxplot Distribusi Nilai Ekspresi",
        ylab = "Expression Value (log2)", cex.axis = 0.6)
legend("topright", legend = levels(gset_sub$group), fill = c("salmon", "turquoise"), cex = 0.8)

# --- 2. DENSITY PLOT ---
library(ggplot2)
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset_sub$group, each = nrow(ex))
)
plot_density <- ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(title = "Distribusi Nilai Ekspresi Gen (Glioblastoma vs Normal)",
       x = "Expression Value (log2)", y = "Density")
print(plot_density)

# --- 3. UMAP PLOT ---
library(umap)
umap_input <- t(ex)
umap_result <- umap(umap_input)
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset_sub$group
)
plot_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP Plot: Glioblastoma vs Normal", x = "UMAP 1", y = "UMAP 2")
print(plot_umap)

# --- PERBAIKAN MUTLAK: Mendefinisikan ulang grup agar R tidak bingung ---
gset_sub$group <- factor(ifelse(grepl("normal", pData(gset_sub)[["source_name_ch1"]], ignore.case=TRUE), "Normal", "Glioblastoma"), levels = c("Glioblastoma", "Normal"))

# --- 1. BOXPLOT WARNA-WARNI (ENGLISH VERSION) ---
group_colors <- ifelse(gset_sub$group == "Glioblastoma", "salmon", "turquoise")
boxplot(ex, col = group_colors, las = 2, outline = FALSE,
        main = "Boxplot of Expression Value Distribution",
        ylab = "Expression Value (log2)", cex.axis = 0.6)
legend("topright", legend = levels(gset_sub$group), fill = c("salmon", "turquoise"), cex = 0.8)

# --- 2. DENSITY PLOT (ENGLISH VERSION) ---
library(ggplot2)
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset_sub$group, each = nrow(ex))
)
plot_density <- ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(title = "Density Plot of Gene Expression (Glioblastoma vs Normal)",
       x = "Expression Value (log2)", y = "Density")
print(plot_density)

# --- 3. UMAP PLOT (ENGLISH VERSION) ---
library(umap)
umap_input <- t(ex)
umap_result <- umap(umap_input)
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset_sub$group
)
plot_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP Plot: Glioblastoma vs Normal", x = "UMAP 1", y = "UMAP 2")
print(plot_umap)