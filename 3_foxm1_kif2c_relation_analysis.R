## Load Packages ---------------------------------------------------------------

set.seed(123)

library(tidyverse)
library(org.Hs.eg.db)
library(patchwork)

## Get TPM From Aging HDFs

aging_133_tpm <- read.csv("input_data\\aging_133_tpm.csv", 
                          header = T, 
                          sep = ";", 
                          row.names = 1) ## Get TPM

metadata <- read.csv("input_data\\sample_metadata_correct_tpm.csv", 
                     header = T, sep = ";") ## Get Metadata

metadata$Sample == colnames(aging_133_tpm)


##Log2(TPM + 1) Tranformation --------------------------------------------------

gene_tpm_log2 <- log2(aging_133_tpm + 1)
ages <- metadata$Age

correlations_log2 <- numeric(nrow(gene_tpm_log2))
pval_log2 <- numeric(nrow(gene_tpm_log2))

## Create loop to obtain Log2(TPM+1) correlation with age:

for (i in 1:nrow(gene_tpm_log2)) {
  counts_log <- as.numeric(gene_tpm_log2[i, ])
  test_log <- cor.test(counts_log, ages, method = "spearman")
  correlations_log2[i] <- test_log$estimate
  pval_log2[i] <- test_log$p.value
} ## This may take a while...

result_df_log2 <- data.frame(
  Gene = row.names(gene_tpm_log2),
  Correlation = correlations_log2,
  P_value = pval_log2
)

result_df_log2 <- result_df_log2 %>% na.omit()

genes_2 <- result_df_log2[, "Gene"]

annots <- AnnotationDbi::select(org.Hs.eg.db, keys=genes_2, 
                                columns="SYMBOL", keytype="ENTREZID")
result_df_log2_2 <- merge(result_df_log2, annots, by.x="Gene", by.y="ENTREZID")

## Now correct for FDR:

result_df_log2_2$fdr <- p.adjust(result_df_log2_2$P_value, method = "BH")

write.csv(result_df_log2_2, "output_data\\excel_results\\aging_log2tpm_spearman_corr.csv")

## Get FOXM1 and KIF2C data:

goi <- c("FOXM1", "KIF2C")

annots <- AnnotationDbi::select(org.Hs.eg.db, keys=goi, 
                                columns="ENTREZID", keytype="SYMBOL")

gene_tpm_log2_goi <- gene_tpm_log2[annots$ENTREZID,]

gene_tpm_log2_goi <- t(gene_tpm_log2_goi) %>% as.data.frame()

gene_tpm_log2_goi$sample <- row.names(gene_tpm_log2_goi)

gene_tpm_log2_goi_2 <- merge(gene_tpm_log2_goi, metadata, by.x = "sample", by.y = "Sample")

colnames(gene_tpm_log2_goi_2)

gene_tpm_log2_goi_2 <- gene_tpm_log2_goi_2 %>% dplyr::rename("FOXM1" = "2305", "KIF2C" = "11004")


gene_tpm_log2_goi_2 <- gene_tpm_log2_goi_2 %>% pivot_longer(names_to = "Gene", values_to = "log2TPM", cols = c(FOXM1, KIF2C))

cols = c(
  "FOXM1" = "gray20",
  "KIF2C" = "darkorange"
)

p1 <- ggplot(gene_tpm_log2_goi_2, aes(x = Age, y = log2TPM, colour = Gene)) + geom_point(size = 2, alpha = 0.7) +
  geom_smooth(se = T, method = "lm") +
  xlim(0, 101) +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 16) +
  labs(x = "Donor Age", y = "Log2(TPM + 1)",
       title = "Gene Expression vs Donor Age") +
  annotate("text", label = "\u03c1 = -0.41 (padjust = 1.21E-5)", x = 0.2, y = 2.2, hjust = 0, color = "gray20", size = 4) +
  annotate("text", label = "\u03c1 = -0.44 (padjust = 1.40E-6)", x = 0.2, y = 1.7, hjust = 0, color = "darkorange", size = 4) 
p1

## Correlation between KIF2C and FOXM1 themselves:

gene_tpm_log2_goi <- gene_tpm_log2[annots$ENTREZID,]

gene_tpm_log2_goi <- t(gene_tpm_log2_goi) %>% as.data.frame()

gene_tpm_log2_goi$sample <- row.names(gene_tpm_log2_goi)

gene_tpm_log2_goi_2 <- merge(gene_tpm_log2_goi, metadata, by.x = "sample", by.y = "Sample")

colnames(gene_tpm_log2_goi_2)

gene_tpm_log2_goi_2 <- gene_tpm_log2_goi_2 %>% dplyr::rename("FOXM1" = "2305", "KIF2C" = "11004")

min(gene_tpm_log2_goi_2$FOXM1)
min(gene_tpm_log2_goi_2$KIF2C)



test_log <- cor.test(gene_tpm_log2_goi_2$FOXM1, gene_tpm_log2_goi_2$KIF2C, method = "spearman")
correl <- test_log$estimate
pval_correl <- test_log$p.value


p2 <- ggplot(gene_tpm_log2_goi_2, aes(x = FOXM1, y = KIF2C)) + geom_point(aes(color = Age), size = 2, alpha = 0.7) +
  scale_color_gradient2(low = "green", mid = "yellow", high = "magenta", limits = c(0, 100), 
                        midpoint = 50) +
  geom_smooth(method = "lm") +
  xlim(2,8) +
  theme_bw(base_size = 16) +
  labs(x = "FOXM1 Log2(TPM + 1)", y = "KIF2C Log2(TPM + 1)",
       title = "KIF2C vs FOXM1 Gene Expression") +
  annotate("text", label = "\u03c1 = 0.93 (p-value = 1.49E-57)", 
           x = 4.2, y = 1.6, hjust = 0, color = "gray20", size = 4)
p2

p3 <- p1 + p2
p3

p3 <- p3 + theme(text = element_text(family = "Arial Unicode MS"))
p3

ggsave(
  "output_data\\plots\\supp_figures\\KIF2C_vs_FOXM1_age_plot.pdf",
  width = 13,
  height = 6,
  unit = "in",
  device = cairo_pdf
)
