library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomicFeatures)
library(dplyr)

## Load DESeq2 Results ---------------------------------------------------------

df_hgps_umk57_vs_dmso <- read.csv(
  "output_data\\excel_results\\deseq2_hgps_umk57_vs_dmso.csv",
  sep = ","
)

df_ds_umk57_vs_dmso <- read.csv(
  "output_data\\excel_results\\deseq2_ds_umk57_vs_dmso.csv",
  sep = ","
) %>% 
  dplyr::mutate(log10padj = -log10(padj)) %>% 
  dplyr::mutate(gene = X)

df_hgps_vs_neo <- read.csv(
  "output_data\\excel_results\\deseq2_hgps_vs_neo.csv",
  sep = ","
)

df_ds_vs_neo <- read.csv(
  "output_data\\excel_results\\deseq2_ds_vs_neo.csv",
  sep = ","
)

all_genes <- rbind(
  df_hgps_umk57_vs_dmso,
  df_ds_umk57_vs_dmso,
  df_ds_vs_neo,
  df_hgps_vs_neo
) %>% 
  dplyr::pull(gene) %>% 
  as.vector() %>% 
  unique() ## Get all expressed genes in all samples

annots <- AnnotationDbi::select(org.Hs.eg.db, keys = all_genes, 
                                columns="ENTREZID", 
                                keytype="SYMBOL") %>% 
  na.omit()


annots <- annots[!duplicated(annots$SYMBOL),]
annots <- annots[!duplicated(annots$ENTREZID),] %>% 
  dplyr::rename(gene = SYMBOL)

df_hgps_umk57_vs_dmso <- df_hgps_umk57_vs_dmso %>% 
  dplyr::left_join(annots, join_by(gene)) %>% na.omit()

df_ds_umk57_vs_dmso <- df_ds_umk57_vs_dmso %>% 
  dplyr::left_join(annots, join_by(gene)) %>% na.omit()

df_hgps_vs_neo <- df_hgps_vs_neo %>% 
  dplyr::left_join(annots, join_by(gene)) %>% na.omit()

df_ds_vs_neo <- df_ds_vs_neo %>% 
  dplyr::left_join(annots, join_by(gene)) %>% na.omit()


## Get Gene Lengths ------------------------------------------------------------

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

genes_gr <- genes(txdb)

gene_length_df <- data.frame(
  ENTREZID = names(genes_gr),
  gene_length = width(genes_gr)
)

gene_length_df$log10gene_length <- log10(gene_length_df$gene_length)

df_hgps_umk57_vs_dmso <- df_hgps_umk57_vs_dmso %>%
  left_join(gene_length_df,
            by = "ENTREZID") %>%
  filter(!is.na(gene_length),
         !is.na(log2FoldChange)) %>% 
  dplyr::mutate(gene_label = case_when(
    gene %in% c("KIF2C",
                "CNTNAP2", "DMD", "ERBB4", "OPCML") ~ gene,
    .default = ""
  ))

df_ds_umk57_vs_dmso <- df_ds_umk57_vs_dmso %>%
  left_join(gene_length_df,
            by = "ENTREZID") %>%
  filter(!is.na(gene_length),
         !is.na(log2FoldChange)) %>% 
  dplyr::mutate(gene_label = case_when(
    gene %in% c("KIF2C",
                "CNTNAP2", "DMD", "ERBB4", "OPCML") ~ gene,
    .default = ""
  ))

df_hgps_vs_neo <- df_hgps_vs_neo %>%
  left_join(gene_length_df,
            by = "ENTREZID") %>%
  filter(!is.na(gene_length),
         !is.na(log2FoldChange)) %>% 
  dplyr::mutate(gene_label = case_when(
    gene %in% c("KIF2C",
                "CNTNAP2", "DMD", "ERBB4", "OPCML") ~ gene,
    .default = ""
  ))

df_ds_vs_neo <- df_ds_vs_neo %>%
  left_join(gene_length_df,
            by = "ENTREZID") %>%
  filter(!is.na(gene_length),
         !is.na(log2FoldChange)) %>% 
  dplyr::mutate(gene_label = case_when(
    gene %in% c("KIF2C",
                "PTPRD", "SGCD", "TENM2", "OPCML") ~ gene,
    .default = ""
  ))


## Spearman Correlation

cor_test_hgps_vs_neo <- cor.test(
  df_hgps_vs_neo$log2FoldChange,
  df_hgps_vs_neo$log10gene_length,
  method = "spearman"
)

cor_test_hgps_vs_neo

spearman_hgps <- cor_test_hgps_vs_neo$estimate
pvalue_hgps <- cor_test_hgps_vs_neo$p.value

cor_test_ds_vs_neo <- cor.test(
  df_ds_vs_neo$log2FoldChange,
  df_ds_vs_neo$log10gene_length,
  method = "spearman"
)

cor_test_ds_vs_neo

spearman_ds <- cor_test_ds_vs_neo$estimate
pvalue_ds <- cor_test_ds_vs_neo$p.value

library(ggpointdensity)

y_limits <- c(10.5, NA)

p1 <- ggplot(df_hgps_vs_neo,
       aes(x = log10gene_length,
           y = log2FoldChange,
           label = gene_label)) +
  geom_pointdensity(adjust = 1, alpha = 0.5) +
  geom_point(data = df_hgps_vs_neo %>% dplyr::filter(gene_label != ""),
             aes(x = log10gene_length, y = log2FoldChange),
             color = "red", size = 3) +
  viridis::scale_color_viridis() +
  ggrepel::geom_text_repel(force_pull   = 0,
                           direction    = "x",
                           angle        = 90,
                           hjust        = 0,
                           ylim = y_limits,
                           segment.size = 0.2,
                           max.iter = 1e4, max.time = 1,
                           color = "black",
                           max.overlaps = Inf) +
  ylim(-11, 14) +
  theme_bw(base_size = 14) +
  annotate("text", x = 2, y = -7, hjust = 0, 
           label = paste0("Rho = ", round(spearman_hgps, 2))) +
  annotate("text", x = 2, y = -8, hjust = 0, 
           label = paste0("pval = ", signif(pvalue_hgps, 3))) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  labs(x = "\nLog10(Gene Length)", y = "LOG2FC\n",
       title = "LOG2FC vs Gene Length - HGPS")

p1 <- ggrastr::rasterise(p1, layers = "Point", dpi = 300)
p1

p2 <- ggplot(df_ds_vs_neo,
       aes(x = log10gene_length,
           y = log2FoldChange,
           label = gene_label)) +
  geom_pointdensity(adjust = 1, alpha = 0.5) +
  geom_point(data = df_ds_vs_neo %>% dplyr::filter(gene_label != ""),
             aes(x = log10gene_length, y = log2FoldChange),
             color = "red", size = 3) +
  viridis::scale_color_viridis() +
  ylim(-11, 14) +
  ggrepel::geom_text_repel(force_pull   = 0,
                           direction    = "x",
                           angle        = 90,
                           hjust        = 0,
                           ylim = y_limits,
                           segment.size = 0.2,
                           max.iter = 1e4, max.time = 1,
                           color = "black",
                           max.overlaps = Inf) + 
  theme_bw(base_size = 14) +
  annotate("text", x = 2, y = -7, hjust = 0, 
           label = paste0("Rho = ", round(spearman_ds, 2))) +
  annotate("text", x = 2, y = -8, hjust = 0, 
           label = paste0("pval = ", signif(pvalue_ds, 3))) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  labs(x = "\nLog10(Gene Length)", y = "LOG2FC\n",
       title = "LOG2FC vs Gene Length - DS")

p2 <- ggrastr::rasterise(p2, layers = "Point", dpi = 300)
p2

library(patchwork)

p1 + p2

ggsave("log2fc_vs_gene_length_disease_vs_neo.pdf",
       path = "output_data\\transcript_length",
       width = 14,
       height = 6)

ggsave("log2fc_vs_gene_length_disease_vs_neo.png",
       path = "output_data\\transcript_length",
       width = 14,
       height = 6)



## Correlation vs UMK57 vs DMSO ------------------------------------------------

## Spearman Correlation

cor_test_hgps_umk57_vs_dmso <- cor.test(
  df_hgps_umk57_vs_dmso$log2FoldChange,
  df_hgps_umk57_vs_dmso$log10gene_length,
  method = "spearman"
)

cor_test_hgps_umk57_vs_dmso

spearman_hgps_umk <- cor_test_hgps_umk57_vs_dmso$estimate
pvalue_hgps_umk <- cor_test_hgps_umk57_vs_dmso$p.value


cor_test_ds_umk57_vs_dmso <- cor.test(
  df_ds_umk57_vs_dmso$log2FoldChange,
  df_ds_umk57_vs_dmso$log10gene_length,
  method = "spearman"
)

cor_test_ds_umk57_vs_dmso

spearman_ds_umk <- cor_test_ds_umk57_vs_dmso$estimate
pvalue_ds_umk <- cor_test_ds_umk57_vs_dmso$p.value


library(ggpointdensity)

y_limits <- c(10.5, NA)

p1 <- ggplot(df_hgps_umk57_vs_dmso,
             aes(x = log10gene_length,
                 y = log2FoldChange,
                 label = gene_label)) +
  geom_pointdensity(adjust = 1, alpha = 0.5) +
  geom_point(data = df_hgps_umk57_vs_dmso %>% dplyr::filter(gene_label != ""),
             aes(x = log10gene_length, y = log2FoldChange),
             color = "red", size = 3) +
  viridis::scale_color_viridis() +
  ggrepel::geom_text_repel(force_pull   = 0,
                           direction    = "x",
                           angle        = 90,
                           hjust        = 0,
                           ylim = y_limits,
                           segment.size = 0.2,
                           max.iter = 1e4, max.time = 1,
                           color = "black",
                           max.overlaps = Inf) +
  ylim(-11, 14) +
  theme_bw(base_size = 14) +
  annotate("text", x = 2, y = -7, hjust = 0, 
           label = paste0("Rho = ", round(spearman_hgps_umk, 2))) +
  annotate("text", x = 2, y = -8, hjust = 0, 
           label = paste0("pval = ", signif(pvalue_hgps_umk, 3))) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  labs(x = "\nLog10(Gene Length)", y = "LOG2FC\n",
       title = "LOG2FC vs Gene Length - HGPS UMK57 vs DMSO")
p1

p1 <- ggrastr::rasterise(p1, layers = "Point", dpi = 300)
p1

p2 <- ggplot(df_ds_umk57_vs_dmso,
             aes(x = log10gene_length,
                 y = log2FoldChange,
                 label = gene_label)) +
  geom_pointdensity(adjust = 1, alpha = 0.5) +
  geom_point(data = df_ds_umk57_vs_dmso %>% dplyr::filter(gene_label != ""),
             aes(x = log10gene_length, y = log2FoldChange),
             color = "red", size = 3) +
  viridis::scale_color_viridis() +
  ylim(-11, 14) +
  ggrepel::geom_text_repel(force_pull   = 0,
                           direction    = "x",
                           angle        = 90,
                           hjust        = 0,
                           ylim = y_limits,
                           segment.size = 0.2,
                           max.iter = 1e4, max.time = 1,
                           color = "black",
                           max.overlaps = Inf) + 
  theme_bw(base_size = 14) +
  annotate("text", x = 2, y = -7, hjust = 0, 
           label = paste0("Rho = ", round(spearman_ds_umk, 2))) +
  annotate("text", x = 2, y = -8, hjust = 0, 
           label = paste0("pval = ", signif(pvalue_ds_umk, 3))) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  labs(x = "\nLog10(Gene Length)", y = "LOG2FC\n",
       title = "LOG2FC vs Gene Length - DS UMK57 vs DMSO")

p2 <- ggrastr::rasterise(p2, layers = "Point", dpi = 300)
p2

library(patchwork)

p1 + p2

ggsave("log2fc_vs_gene_length_umk57_vs_dmso.pdf",
       path = "output_data\\transcript_length",
       width = 14,
       height = 6)

ggsave("log2fc_vs_gene_length_umk57_vs_dmso.png",
       path = "output_data\\transcript_length",
       width = 14,
       height = 6)

