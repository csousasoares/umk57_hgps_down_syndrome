## Load Packages and Set Seed --------------------------------------------------

set.seed(123)

library(DESeq2)
library(org.Hs.eg.db)
library(tidyverse)
library(stringr)
library(ggrepel)
library(ggpubr)
library(ggalt)
library(limma)
library(circlize)
library(ComplexHeatmap)
library(colorRamp2)
library(ggtree)
library(patchwork)

## Load Counts and Metadata ----------------------------------------------------

## Load correct sample info and counts for DS, HGPS and Healthy
## samples obtained from IonTorrent Suite:

counts_matrix_correct <- read.csv("input_data\\counts_matrix.csv",
                                  row.names = 1)

sample_info_correct <- read.csv("input_data\\sample_info.csv", 
                                row.names = 1)



## DESeq2 HGPS vs Healthy ------------------------------------------------------

sample_info_correct_hgps_neo <- sample_info_correct %>% 
  dplyr:: filter(Disease %in% c("Healthy", "HGPS")) %>% 
  dplyr::filter(Treatment == "DMSO")
sample_info_correct_hgps_neo

counts_matrix_correct_hgps_neo <- counts_matrix_correct[,row.names(sample_info_correct_hgps_neo)]

colnames(counts_matrix_correct_hgps_neo) == rownames(sample_info_correct_hgps_neo)
## Should be TRUE

dds_hgps_neo <- DESeqDataSetFromMatrix(
  countData = counts_matrix_correct_hgps_neo,
  colData = sample_info_correct_hgps_neo,
  design = ~ Disease)

## Design wont matter at this stage though

smallestGroupSize <- nrow(sample_info_correct_hgps_neo)/2 
## Half of samples
keep <- rowSums(counts(dds_hgps_neo) >=5) >= smallestGroupSize 
##Recomended in vignette
dds_hgps_neo <- dds_hgps_neo[keep,]

dds_hgps_neo

dds_hgps_neo <- DESeq(dds_hgps_neo)
res_hgps_neo <- results(dds_hgps_neo, contrast = 
                          c("Disease", "HGPS", "Healthy"), 
                        cooksCutoff=FALSE, independentFiltering=FALSE) 


df_hgps_neo <- as.data.frame(res_hgps_neo)

## GSEA HGPS vs Healthy Hallmark -----------------------------------------------

head(df_hgps_neo)

df_hgps_neo <- df_hgps_neo %>% 
  dplyr::arrange(desc(log2FoldChange))
head(df_hgps_neo)

df_hgps_neo$gene <- rownames(df_hgps_neo)

write.csv(df_hgps_neo, "output_data\\excel_results\\deseq2_hgps_vs_neo.csv")

hgps_neo_ranked <- df_hgps_neo$log2FoldChange
names(hgps_neo_ranked) <- df_hgps_neo$gene
head(hgps_neo_ranked, 10)


hall <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hall_t2g <- hall %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

hgps_neo_hall <- clusterProfiler::GSEA(
  hgps_neo_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = hall_t2g,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

df_hgps_neo_hall <- as.data.frame(hgps_neo_hall)

df_hgps_neo_hall_005 <- df_hgps_neo_hall %>% 
  dplyr::filter(p.adjust < 0.05)

df_hgps_neo_hall_005$Description <- gsub(
  "HALLMARK_", "", df_hgps_neo_hall_005$Description
)

df_hgps_neo_hall_005 <- df_hgps_neo_hall_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

cols <- c(
  "Up" = "red",
  "Down" = "blue"
)

ggplot(df_hgps_neo_hall_005,
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", 
                        name = "padjust") +
  theme_bw(base_size = 13) +
  xlim(-4, 3) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "HGPS vs NEO",
       subtitle = "Hallmark Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")


## DESeq2 DS vs Healthy --------------------------------------------------------

sample_info_correct_ds_neo <- sample_info_correct %>% 
  dplyr::filter(Treatment == "DMSO") %>% 
  dplyr::filter(Disease %in% c("Healthy", "DS"))

counts_matrix_correct_ds_neo <- counts_matrix_correct[,row.names(sample_info_correct_ds_neo)]


dds_ds_neo <- DESeqDataSetFromMatrix(
  countData = counts_matrix_correct_ds_neo,
  colData = sample_info_correct_ds_neo,
  design = ~ Disease) 

## Design wont matter at this stage though

smallestGroupSize <- nrow(sample_info_correct_ds_neo)/2 
## Half of samples, can be changed later
keep <- rowSums(counts(dds_ds_neo) >=5) >= smallestGroupSize 
##Recomended in vignette
dds_ds_neo <- dds_ds_neo[keep,]

dds_ds_neo

dds_ds_neo <- DESeq(dds_ds_neo)
res_ds_neo <- results(dds_ds_neo, contrast = 
                        c("Disease", "DS", "Healthy"), 
                      cooksCutoff=FALSE, independentFiltering=FALSE) 


df_ds_neo <- as.data.frame(res_ds_neo)

## GSEA DS vs Healthy Hallmark -------------------------------------------------

head(df_ds_neo)

df_ds_neo$gene <- rownames(df_ds_neo)

df_ds_neo <- df_ds_neo %>% 
  dplyr::arrange(desc(log2FoldChange))
head(df_ds_neo)

write.csv(df_ds_neo, "output_data\\excel_results\\deseq2_ds_vs_neo.csv")

ds_neo_ranked <- df_ds_neo$log2FoldChange
names(ds_neo_ranked) <- df_ds_neo$gene
head(ds_neo_ranked, 10)


hall <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hall_t2g <- hall %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

ds_neo_hall <- clusterProfiler::GSEA(
  ds_neo_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = hall_t2g,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

df_ds_neo_hall <- as.data.frame(ds_neo_hall)

df_ds_neo_hall_005 <- df_ds_neo_hall %>% 
  dplyr::filter(p.adjust < 0.05)

df_ds_neo_hall_005$Description <- gsub(
  "HALLMARK_", "", df_ds_neo_hall_005$Description
)

df_ds_neo_hall_005 <- df_ds_neo_hall_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

cols <- c(
  "Up" = "red",
  "Down" = "blue"
)

ggplot(df_ds_neo_hall_005,
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", 
                        name = "padjust") +
  theme_bw(base_size = 13) +
  xlim(-4,3) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "DS vs NEO",
       subtitle = "Hallmark Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")


## Combined Dotplot Hallmark ---------------------------------------------------

df_ds_neo_hall[1:5,2:5]
df_hgps_neo_hall[1:5,2:5]
df_hgps_umk57_hall <-  read.csv("output_data\\excel_results\\hgps_umk57_dmso_hallmarks.csv",
                                row.names = 1)
df_Hall_ds_of_interest <- read.csv("output_data\\excel_results\\ds_umk57_dmso_hallmarks.csv",
                                   row.names = 1)

df_ds_neo_hall$comparison <- "DS vs Healthy"
df_hgps_neo_hall$comparison <- "HGPS vs Healthy"
df_hgps_umk57_hall$comparison <- "HGPS UMK57 vs DMSO"
df_Hall_ds_of_interest$comparison <- "DS UMK57 vs DMSO"

colnames(df_ds_neo_hall) == colnames(df_Hall_ds_of_interest)

dodgeplot_ds <- rbind(df_ds_neo_hall,
                      df_Hall_ds_of_interest)

dodgeplot_ds$Description <- gsub(
  "HALLMARK_",
  "",
  dodgeplot_ds$Description
)

dodgeplot_ds <- dodgeplot_ds %>% 
  dplyr::mutate(NES_NEO = case_when(
    comparison == "DS vs Healthy" ~ NES,
    .default = 0
  ))

dodgeplot_ds <- dodgeplot_ds %>% 
  dplyr::arrange(desc(NES_NEO))

dodgeplot_ds_005 <- dodgeplot_ds %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::arrange(desc(NES_NEO))

unique(dodgeplot_ds_005$Description)

dodgeplot_hgps <- rbind(df_hgps_neo_hall,
                        df_hgps_umk57_hall)

dodgeplot_hgps$Description <- gsub(
  "HALLMARK_",
  "",
  dodgeplot_hgps$Description
)

dodgeplot_hgps <- dodgeplot_hgps %>% 
  dplyr::mutate(NES_NEO = case_when(
    comparison == "HGPS vs Healthy" ~ NES,
    .default = 0
  ))

dodgeplot_hgps <- dodgeplot_hgps %>% 
  dplyr::arrange(desc(NES_NEO))

dodgeplot_hgps_005 <- dodgeplot_hgps %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::arrange(desc(NES_NEO))

unique(dodgeplot_hgps_005$Description)


dodgeplot_hall_combined <- rbind(
  dodgeplot_ds_005,
  dodgeplot_hgps_005
)

dodgeplot_hall_combined$comparison <- factor(
  dodgeplot_hall_combined$comparison,
  levels = c("HGPS vs Healthy",
             "HGPS UMK57 vs DMSO",
             "DS vs Healthy",
             "DS UMK57 vs DMSO")
)

hall_factors <- dodgeplot_hall_combined %>% 
  filter(NES_NEO != 0) %>% 
  arrange(desc(NES_NEO)) %>% 
  pull(Description) %>% unique() %>% as.vector()

hall_factors_2 <- dodgeplot_hall_combined %>% 
  pull(Description) %>% unique() %>% as.vector()

hall_factors_2 <- setdiff(hall_factors_2, hall_factors)

hall_factors_final <- c(hall_factors_2, hall_factors)

dodgeplot_hall_combined$Description <- factor(
  dodgeplot_hall_combined$Description,
  levels = hall_factors_final
)

ggplot(dodgeplot_hall_combined, aes(
  x = comparison, 
  y = Description, 
  group = comparison)) +
  geom_point(shape = 21, aes(fill = NES, size = p.adjust)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       breaks = c(-3:2)) +
  scale_size_continuous(range = c(7, 2.5), transform = "log10",
                        name = "padjust") +
  theme_bw(base_size = 14) +
  labs(y = "Gene Sets\n", x = "",
       title = "Hallmark Gene Sets") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_line(size = 0.25)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(fill = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.25,
    ticks.linewidth = 0.25,
    ticks.colour = "black")) +
  geom_vline(xintercept = 2.5, linetype = "dotted")

ggsave("hall_dotplot_all.pdf",
       path = "output_data\\plots\\supp_figures",
       width = 7,
       height = 10)


## GSEA HGPS vs Healthy DNA Repair ---------------------------------------------

mm_BP_sets <- msigdbr::msigdbr(
  species = "Homo sapiens")
mm_BP_sets

msigdbr_t2g = mm_BP_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% 
  as.data.frame()
msigdbr_t2g

misgdbr_t2g_filtered <- filter(msigdbr_t2g, gs_name %in% c("GOBP_BASE_EXCISION_REPAIR",
                                                           "GOBP_DNA_REPAIR",
                                                           "GOBP_DNA_SYNTHESIS_INVOLVED_IN_DNA_REPAIR",
                                                           "GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
                                                           "GOBP_INTERSTRAND_CROSS_LINK_REPAIR",
                                                           "GOBP_MISMATCH_REPAIR",
                                                           "GOBP_NUCLEOTIDE_EXCISION_REPAIR",
                                                           "GOBP_POSITIVE_REGULATION_OF_DNA_REPAIR",
                                                           "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
                                                           "GOBP_POSTREPLICATION_REPAIR",
                                                           "GOBP_RECOMBINATIONAL_REPAIR",
                                                           "GOBP_REGULATION_OF_DNA_REPAIR",
                                                           "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
                                                           "HALLMARK_DNA_REPAIR",
                                                           "KAUFFMANN_DNA_REPAIR_GENES",
                                                           "KEGG_BASE_EXCISION_REPAIR",
                                                           "KEGG_NUCLEOTIDE_EXCISION_REPAIR",
                                                           "REACTOME_BASE_EXCISION_REPAIR",
                                                           "REACTOME_BASE_EXCISION_REPAIR_AP_SITE_FORMATION",
                                                           "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "REACTOME_DNA_REPAIR",
                                                           "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
                                                           "REACTOME_HOMOLOGY_DIRECTED_REPAIR",
                                                           "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
                                                           "REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS",
                                                           "REACTOME_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR_TC_NER",
                                                           "WP_BASE_EXCISION_REPAIR",
                                                           "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
                                                           "WP_NUCLEOTIDE_EXCISION_REPAIR"))
misgdbr_t2g_filtered


hgps__vs_neo_dna_repair <- clusterProfiler::GSEA(
  hgps_neo_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = misgdbr_t2g_filtered,
  verbose = T,
  seed = T,
  by = "fgsea",
  nPerm = 20000 ## Since it shows warning, to make sure all terms are included
)

df_hgps__vs_neo_dna_repair <- as.data.frame(hgps__vs_neo_dna_repair)

df_hgps__vs_neo_dna_repair_005 <- df_hgps__vs_neo_dna_repair %>% 
  dplyr::filter(p.adjust < 0.05)

df_hgps__vs_neo_dna_repair_005$Description <- gsub("_", " ", df_hgps__vs_neo_dna_repair_005$Description)

df_hgps__vs_neo_dna_repair_005$Description <- stringr::str_wrap(df_hgps__vs_neo_dna_repair_005$Description,
                                                                width = 30)

df_hgps__vs_neo_dna_repair_005$Description

df_hgps__vs_neo_dna_repair_005 <- df_hgps__vs_neo_dna_repair_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

cols = c(
  "Up" = "red",
  "Down" = "blue"
)

ggplot(df_hgps__vs_neo_dna_repair_005,
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", name = "padjust") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(-3,3) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "HGPS vs NEO",
       subtitle = "DNA Repair Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")


## GSEA DS vs Healthy DNA Repair -----------------------------------------------

ds__vs_neo_dna_repair <- clusterProfiler::GSEA(
  ds_neo_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = misgdbr_t2g_filtered,
  verbose = T,
  seed = T,
  by = "fgsea",
)


df_ds__vs_neo_dna_repair <- as.data.frame(ds__vs_neo_dna_repair)

df_ds__vs_neo_dna_repair_005 <- df_ds__vs_neo_dna_repair %>% 
  dplyr::filter(p.adjust < 0.05)

df_ds__vs_neo_dna_repair_005$Description <- gsub("_", " ", df_ds__vs_neo_dna_repair_005$Description)

df_ds__vs_neo_dna_repair_005$Description <- stringr::str_wrap(df_ds__vs_neo_dna_repair_005$Description,
                                                              width = 30)

df_ds__vs_neo_dna_repair_005$Description

df_ds__vs_neo_dna_repair_005 <- df_ds__vs_neo_dna_repair_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

ggplot(df_ds__vs_neo_dna_repair_005,
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", name = "padjust") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(-4,3) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "DS vs NEO",
       subtitle = "DNA Repair Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")



## Divergent Dotplots DNA Repair -----------------------------------------------

new_dna_repair_terms <- data.frame(
  ID = c("GOBP_BASE_EXCISION_REPAIR",
         "GOBP_DNA_REPAIR",
         "GOBP_DNA_SYNTHESIS_INVOLVED_IN_DNA_REPAIR",
         "GOBP_DOUBLE_STRAND_BREAK_REPAIR",
         "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
         "GOBP_INTERSTRAND_CROSS_LINK_REPAIR",
         "GOBP_MISMATCH_REPAIR",
         "GOBP_NUCLEOTIDE_EXCISION_REPAIR",
         "GOBP_POSITIVE_REGULATION_OF_DNA_REPAIR",
         "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR",
         "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
         "GOBP_POSTREPLICATION_REPAIR",
         "GOBP_RECOMBINATIONAL_REPAIR",
         "GOBP_REGULATION_OF_DNA_REPAIR",
         "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR",
         "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
         "HALLMARK_DNA_REPAIR",
         "KAUFFMANN_DNA_REPAIR_GENES",
         "KEGG_BASE_EXCISION_REPAIR",
         "KEGG_NUCLEOTIDE_EXCISION_REPAIR",
         "REACTOME_BASE_EXCISION_REPAIR",
         "REACTOME_BASE_EXCISION_REPAIR_AP_SITE_FORMATION",
         "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
         "REACTOME_DNA_REPAIR",
         "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
         "REACTOME_HOMOLOGY_DIRECTED_REPAIR",
         "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
         "REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS",
         "REACTOME_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR_TC_NER",
         "WP_BASE_EXCISION_REPAIR",
         "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
         "WP_NUCLEOTIDE_EXCISION_REPAIR"),
  "new" = c(1:32))


new_dna_repair_terms$new <- c(
  "(BP) Base Excision Repair",
  "(BP) DNA Repair",
  "(BP) DNA Synthesis Involved in DNA Repair",
  "(BP) DSB Repair",
  "(BP) DSB Repair via NHEJ",
  "(BP) Interstrand Crosslink Repair",
  "(BP) Mismatch Repair",
  "(BP) Nucleotide Excision Repair",
  "(BP) Positive Regulation of DNA Repair",
  "(BP)  Positive Regulation of DSB Repair",
  "(BP) Positive Regulation of DSB Repair via HR",
  "(BP) Postreplication Repair",
  "(BP) Recombinational Repair",
  "(BP) Regulation of DNA Repair",
  "(BP) Regulation of DSB Repair",
  "(BP) Regulation of DSB Repair via HR",
  "(H) DNA Repair",
  "Kauffman DNA Repair Genes",
  "(KEGG) Base Excision Repair",
  "(KEGG) Nucleotide Excision Repair (NER)",
  "(R) Base Excision Repair (BER)",
  "(R) BER - AP Site Formation",
  "(R) DNA DSB Repair",
  "(R) DNA Repair",
  "(R) GG-NER",
  "(R) Homology-Directed Repair",
  "(R) Nucleotide Excision Repair (NER)",
  "(R) Sumoylation of DNA Repair Proteins",
  "(R) TC-NER",
  "(WP) Base Excision Repair (BER)",
  "(WP) DNA Repair Pathways Full Network",
  "(WP) Nucleotide Excision Repair (NER)"
)


df_t21_dna_repair_005 <- readRDS("output_data\\df_t21_dna_repair_005.rds")
df_hgps_dna_repair_005 <- readRDS("output_data\\df_hgps_dna_repair_005.rds")
df_hgps__vs_neo_dna_repair_005
df_ds__vs_neo_dna_repair_005


df_hgps__vs_neo_dna_repair_005$comparison <- "HGPS vs Healthy Control"
df_ds__vs_neo_dna_repair_005$comparison <- "DS vs Healthy Control"

ds_dna_repair_divergent <- rbind(
  df_t21_dna_repair_005,
  df_ds__vs_neo_dna_repair_005
)

ds_dna_repair_divergent$comparison <- factor(
  ds_dna_repair_divergent$comparison,
  levels = c("DS vs Healthy Control",
             "DS UMK57 vs DMSO")
)

ds_dna_repair_divergent$Description <- str_wrap(
  ds_dna_repair_divergent$Description,
  width = 900) ## Make sure no wraping is made now.

df_ds__vs_neo_dna_repair_005$Description <- str_wrap(
  df_ds__vs_neo_dna_repair_005$Description,
  width = 900) ## Make sure no wraping is made now.

ds_dna_repair_divergent <- ds_dna_repair_divergent %>% 
  left_join(new_dna_repair_terms, join_by(ID))

ds_order <- ds_dna_repair_divergent %>% 
  dplyr::filter(comparison == "DS vs Healthy Control") %>% 
  dplyr::arrange(desc(NES)) %>% 
  dplyr::pull(new)

ds_dna_repair_divergent$new <- factor(
  ds_dna_repair_divergent$new,
  levels = ds_order
)

ds_neo_dna_rep_plot <- ggplot(ds_dna_repair_divergent, aes(
  x = comparison, 
  y = new, 
  group = comparison)) +
  geom_point(shape = 21, aes(fill = NES, size = p.adjust)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  scale_size_continuous(range = c(6, 2.5), transform = "log10",
                        name = "padjust") +
  theme_bw(base_size = 14) +
  labs(y = "Gene Sets\n", x = "",
       title = "DS - DNA Repair Gene Sets") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        panel.grid.major = element_line(size = 0.25)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(fill = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.25,
    ticks.linewidth = 0.25,
    ticks.colour = "black"))

ds_neo_dna_rep_plot

hgps_dna_repair_divergent <- rbind(
  df_hgps_dna_repair_005,
  df_hgps__vs_neo_dna_repair_005
)

hgps_dna_repair_divergent$comparison <- factor(
  hgps_dna_repair_divergent$comparison,
  levels = c("HGPS vs Healthy Control",
             "HGPS UMK57 vs DMSO")
)

hgps_dna_repair_divergent$Description <- str_wrap(
  hgps_dna_repair_divergent$Description,
  width = 900) ## Make sure no wraping is made now.

df_hgps__vs_neo_dna_repair_005$Description <- str_wrap(
  df_hgps__vs_neo_dna_repair_005$Description,
  width = 900) ## Make sure no wraping is made now.

hgps_dna_repair_divergent <- hgps_dna_repair_divergent %>% 
  left_join(new_dna_repair_terms, join_by(ID))

hgps_order <- hgps_dna_repair_divergent %>% 
  dplyr::filter(comparison == "HGPS vs Healthy Control") %>% 
  dplyr::arrange(desc(NES)) %>% 
  dplyr::pull(new)

hgps_dna_repair_divergent$new <- factor(
  hgps_dna_repair_divergent$new,
  levels = hgps_order
)

hgps_neo_dna_rep_plot <- ggplot(hgps_dna_repair_divergent, aes(
  x = comparison, 
  y = new, 
  group = comparison)) +
  geom_point(shape = 21, aes(fill = NES, size = p.adjust)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  scale_size_continuous(range = c(6, 2.5), transform = "log10",
                        name = "padjust") +
  theme_bw(base_size = 14) +
  labs(y = "Gene Sets\n", x = "",
       title = "HGPS - DNA Repair Gene Sets") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        panel.grid.major = element_line(size = 0.25)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(fill = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.25,
    ticks.linewidth = 0.25,
    ticks.colour = "black"))

hgps_neo_dna_rep_plot

hgps_neo_dna_rep_plot + ds_neo_dna_rep_plot

ggsave("dotplots_dna_repair_divergent.pdf",
       path = "output_data\\plots\\plots_for_publication",
       width = 12,
       height = 8)


all_dna_repair_divergent <- rbind(
  hgps_dna_repair_divergent,
  ds_dna_repair_divergent
)

all_dna_repair_divergent$comparison <- factor(
  all_dna_repair_divergent$comparison,
  levels = c("HGPS vs Healthy Control",
             "HGPS UMK57 vs DMSO",
             "DS vs Healthy Control",
             "DS UMK57 vs DMSO")
)

ggplot(all_dna_repair_divergent, aes(
  x = comparison, 
  y = forcats::fct_reorder(new, setSize), 
  group = comparison)) +
  geom_point(shape = 21, aes(fill = NES, size = p.adjust)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  scale_size_continuous(range = c(6, 2.5), transform = "log10",
                        name = "padjust") +
  theme_bw(base_size = 14) +
  labs(y = "Gene Sets\n", x = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        panel.grid.major = element_line(size = 0.25)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(fill = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.25,
    ticks.linewidth = 0.25,
    ticks.colour = "black"))

ggsave("dna_repair_dotplot_all.pdf",
       path = "output_data\\plots\\supp_figures",
       width = 9,
       height = 10)


## GSEA DNA Repair Running Plots -----------------------------------------------

## DSB Repair Only:

library(enrichplot)

trace("gseaplot2", edit = T)

hgps_dna_repair <- readRDS("output_data\\hgps_dna_repair.rds")
t21_dna_repair <- readRDS("output_data\\t21_dna_repair.rds")

df_hgps_dna_repair <- as.data.frame(hgps_dna_repair)
df_t21_dna_repair <- as.data.frame(t21_dna_repair)


dsb_line <- match("GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                  df_hgps_dna_repair$Description)

dna_rep_cols <- c(
  "HGPS UMK57 vs DMSO" = "#834746",
  "DS UMK57 vs DMSO" = "lightskyblue"
)

p2 <- enrichplot::gseaplot2(hgps_dna_repair, geneSetID = dsb_line,
                            title = "DSB Repair HGPS UMK57 vs DMSO",
                            base_size = 15,
                            color = "#834746",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p2


df_hgps__vs_neo_dna_repair

dsb_line <- match("GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                  df_hgps__vs_neo_dna_repair$Description)

p1 <- enrichplot::gseaplot2(hgps__vs_neo_dna_repair, geneSetID = dsb_line,
                            title = "DSB Repair HGPS vs NEO",
                            base_size = 15,
                            color = "#834746",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p1


df_ds__vs_neo_dna_repair

dsb_line <- match("GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                  df_ds__vs_neo_dna_repair$Description)

p3 <- enrichplot::gseaplot2(ds__vs_neo_dna_repair, geneSetID = dsb_line,
                            title = "DSB Repair DS vs NEO",
                            base_size = 15,
                            color = "lightskyblue",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p3


df_t21_dna_repair

dsb_line <- match("GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                  df_t21_dna_repair$Description)

p4 <- enrichplot::gseaplot2(t21_dna_repair, geneSetID = dsb_line,
                            title = "DSB Repair DS UMK57 vs DMSO",
                            base_size = 15,
                            color = "lightskyblue",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p4


plot_list(p1, p2, p3, p4, nrow = 1)

ggsave("gseaplot2_gobp_dsb_repair_all.pdf",
       path = "output_data\\plots\\supp_figures/",
       width = 26,
       height = 5)


## GSEA SenMayo ----------------------------------------------------------------

## Ranked Lists of Interest

ds_of_interest_ranked <- readRDS("output_data\\ds_of_interest_ranked.rds")
hgps_umk57_ranked <- readRDS("output_data\\hgps_umk57_ranked.rds")

head(ds_neo_ranked)
head(hgps_neo_ranked)
head(ds_of_interest_ranked)
head(hgps_umk57_ranked)

senmayo <- read.table("input_data\\senmayo.txt",
                      header = T) %>% 
  dplyr::filter(gene_set == "SenMayo")

## DS vs NEO

ds_neo_senmayo <- clusterProfiler::GSEA(
  ds_neo_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = senmayo,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

## HGPS vs NEO

hgps_neo_senmayo <- clusterProfiler::GSEA(
  hgps_neo_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = senmayo,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)


## HGPS UMK57 vs DMSO

hgps_umk57_dmso_senmayo <- clusterProfiler::GSEA(
  hgps_umk57_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = senmayo,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

## DS UMK57 vs DMSO

ds_umk57_dmso_senmayo <- clusterProfiler::GSEA(
  ds_of_interest_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = senmayo,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

## GSEA SenMayo Running Plots --------------------------------------------------


p2 <- enrichplot::gseaplot2(hgps_umk57_dmso_senmayo, geneSetID = 1,
                            title = "SenMayo HGPS UMK57 vs DMSO",
                            base_size = 15,
                            color = "#834746",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p2


p1 <- enrichplot::gseaplot2(hgps_neo_senmayo, geneSetID = 1,
                            title = "SenMayo HGPS vs NEO",
                            base_size = 15,
                            color = "#834746",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p1

p3 <- enrichplot::gseaplot2(ds_neo_senmayo, geneSetID = 1,
                            title = "SenMayo DS vs NEO",
                            base_size = 15,
                            color = "lightskyblue",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p3


p4 <- enrichplot::gseaplot2(ds_umk57_dmso_senmayo, geneSetID = 1,
                            title = "SenMayo DS UMK57 vs DMSO",
                            base_size = 15,
                            color = "lightskyblue",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p4


plot_list(p1, p2, p3, p4, nrow = 1)

ggsave("gseaplot2_senmayo_all.pdf",
       path = "output_data\\plots\\supp_figures",
       width = 26,
       height = 5)

## SenMayo Heatmap -------------------------------------------------------------



## Heatmap FOXM1, KIF2C, etc ---------------------------------------------------

## For HGPS:

sample_info_hgps_umk57_neo <- sample_info_correct %>% 
  dplyr::filter(Disease %in% c("HGPS", "Healthy"))
  
counts_matrix_neo_hgps_umk57 <- counts_matrix_correct[,row.names(
  sample_info_hgps_umk57_neo)]

colnames(counts_matrix_neo_hgps_umk57) == row.names(sample_info_hgps_umk57_neo)

dds_neo_hgps_umk57 <- DESeqDataSetFromMatrix(
  countData = counts_matrix_neo_hgps_umk57,
  colData = sample_info_hgps_umk57_neo,
  design = ~ Individual) ## Makes no difference for PCA

colnames(dds_neo_hgps_umk57)

## Design wont matter at this stage though

smallestGroupSize <- 5 ## Half of samples, can be changed later
keep <- rowSums(counts(dds_neo_hgps_umk57) >=5) >= smallestGroupSize 
##Recomended in vignette
dds_neo_hgps_umk57 <- dds_neo_hgps_umk57[keep,]

dds_neo_hgps_umk57

dds_neo_hgps_umk57 <- DESeq(dds_neo_hgps_umk57)

dds_neo_hgps_umk57

sample_info_hgps_umk57_neo

vsd <- vst(dds_neo_hgps_umk57, blind=TRUE)

vsd_df <- assay(vsd)

sample_info_hgps_umk57_neo <- sample_info_hgps_umk57_neo %>% 
  dplyr::arrange(desc(Individual), Treatment)

vsd_df <- vsd_df[,row.names(sample_info_hgps_umk57_neo)]

vsd_zscore <- t(scale(t(vsd_df)))

goi <- c("FOXM1", "KIF2C", "LMNB1", "IL6", "CDKN1A")

vsd_zscore_goi <- vsd_zscore %>% as.data.frame()
vsd_zscore_goi <- vsd_zscore_goi[goi,]
vsd_zscore_goi <- na.omit(vsd_zscore_goi)

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-2, 2)) ## Necessary to assign values to colors



ha1 = HeatmapAnnotation(Individual = sample_info_hgps_umk57_neo$Individual,
                        Treatment = sample_info_hgps_umk57_neo$Treatment,
                        col = list(
                          Individual = c(
                            "HGPS 169" = "green",
                            "NEO" = "gray"
                          ),
                          Treatment = c(
                            "DMSO" = "gray",
                            "UMK57" = "cyan"
                          )
                        ),
                        
                        gp = gpar(col = "black"))



group_split <- sample_info_hgps_umk57_neo$Disease 

group_split <- factor(group_split,
                      levels = c("Healthy", "HGPS"))


ht <- Heatmap(
  vsd_zscore_goi,
  cluster_columns = F,
  show_column_names = F,
  column_split = group_split,
  top_annotation = ha1,
  col = col_fun,
  row_km = 2,
  border_gp = gpar(col = "black", lty = 1),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2"),
    title = "Row Z-Score",
    border = "black",
    title_position = "leftcenter-rot"))
ht

pdf("output_data\\plots\\supp_figures\\heatmap_healthy_hgps_healthy_foxm1.pdf",
    width = 6,
    height = 4)
draw(ht)
dev.off()


## For DS:

sample_info_neo_good <- sample_info_correct %>% 
  dplyr::filter(Individual == "NEO" & Treatment == "DMSO") %>% 
  dplyr::filter(batch != "07082019") %>% 
  dplyr::select(1:3)

sample_info_ds_umk57_neo <- sample_info_correct %>% 
  dplyr::filter(Disease %in% c("DS", "Healthy"))

counts_matrix_neo_ds_umk57 <- counts_matrix_correct[,row.names(
  sample_info_ds_umk57_neo)]

colnames(counts_matrix_neo_ds_umk57) == row.names(sample_info_ds_umk57_neo)

dds_neo_ds_umk57 <- DESeqDataSetFromMatrix(
  countData = counts_matrix_neo_ds_umk57,
  colData = sample_info_ds_umk57_neo,
  design = ~ Individual) ## Makes no difference for PCA

colnames(dds_neo_ds_umk57)

## Design wont matter at this stage though

smallestGroupSize <- 5 ## Half of samples, can be changed later
keep <- rowSums(counts(dds_neo_ds_umk57) >=5) >= smallestGroupSize 
##Recomended in vignette
dds_neo_ds_umk57 <- dds_neo_ds_umk57[keep,]

dds_neo_ds_umk57

dds_neo_ds_umk57 <- DESeq(dds_neo_ds_umk57)

dds_neo_ds_umk57

sample_info_ds_umk57_neo

vsd <- vst(dds_neo_ds_umk57, blind=TRUE)

vsd_df <- assay(vsd)

sample_info_ds_umk57_neo <- sample_info_ds_umk57_neo %>% 
  dplyr::arrange(Treatment, desc(Individual))

vsd_df <- vsd_df[,row.names(sample_info_ds_umk57_neo)]

vsd_zscore <- t(scale(t(vsd_df)))

goi <- c("FOXM1", "KIF2C", "LMNB1", "IL6", "CDKN1A")

vsd_zscore_goi <- vsd_zscore %>% as.data.frame()
vsd_zscore_goi <- vsd_zscore_goi[goi,]
vsd_zscore_goi <- na.omit(vsd_zscore_goi)

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-2, 2)) ## Necessary to assign values to colors



ha1 = HeatmapAnnotation(Individual = sample_info_ds_umk57_neo$Individual,
                        Treatment = sample_info_ds_umk57_neo$Treatment,
                        col = list(
                          Individual = c(
                            "14y_DS" = "blue",
                            "5y_DS" = "magenta",
                            "NEO" = "gray"
                          ),
                          Treatment = c(
                            "DMSO" = "gray",
                            "UMK57" = "cyan"
                          )
                        ),
                        
                        gp = gpar(col = "black"))



group_split <- sample_info_ds_umk57_neo$Disease 

group_split <- factor(group_split,
                      levels = c("Healthy", "DS"))


ht_2 <- Heatmap(
  vsd_zscore_goi,
  cluster_columns = F,
  show_column_names = F,
  column_split = group_split,
  top_annotation = ha1,
  col = col_fun,
  row_km = 2,
  border_gp = gpar(col = "black", lty = 1),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2"),
    title = "Row Z-Score",
    border = "black",
    title_position = "leftcenter-rot"))
ht_2

pdf("output_data\\plots\\supp_figures\\heatmap_healthy_ds_healthy_foxm1.pdf",
    width = 6,
    height = 4)
draw(ht)
dev.off()


ht_3 <- ht + ht_2
ht_3

## Export Manually 12.5x4


## SenMayo DEGs in NEO + HGPS Full Overlap -------------------------------------

vsd <- vst(dds_neo_hgps_umk57, blind=TRUE) ## Blind to ignore conditions

degs_hgps_vs_neo <- df_hgps_neo %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::pull(gene)

degs_hgps <- read.csv("output_data\\excel_results\\deseq2_hgps_umk57_vs_dmso.csv") %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::pull(gene)

senmayo_genes <- senmayo$gene_symbol

senamyo_degs_hgps_vs_neo <- intersect(
  degs_hgps_vs_neo,
  senmayo_genes
)

senamyo_degs_hgps_vs_neo <- intersect(
  senamyo_degs_hgps_vs_neo,
  degs_hgps
)
## Create heatmap with genes of interest:

vsd_df <- assay(vsd)
vsd_zscore <- scale(t(vsd_df))

vsd_zscore_goi <- vsd_zscore %>% as.data.frame()
vsd_zscore_goi <- vsd_zscore_goi[rownames(hgps_correct_order),senamyo_degs_hgps_vs_neo]
vsd_zscore_goi <- na.omit(vsd_zscore_goi)


col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-2, 2)) ## Necessary to assign values to colors

ha1 = rowAnnotation(Individual = hgps_correct_order$Individual,
                    Treatment = hgps_correct_order$Treatment,
                    col = list(
                      Individual = c(
                        "NEO" = "gray",
                        "HGPS 169" = "green"
                      ),
                      Treatment = c(
                        "DMSO" = "gray",
                        "UMK57" = "cyan"
                      )
                    ))



group_split <- hgps_correct_order$Individual      

group_split <- factor(group_split, levels = c("NEO", "HGPS 169"))


ht_2 <- Heatmap(
  vsd_zscore_goi,
  cluster_columns = T,
  cluster_rows = F,
  show_row_names = F,
  row_split = group_split,
  column_km = 2,
  left_annotation = ha1,
  col = col_fun,
  border_gp = gpar(col = "black", lty = 1),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 90,
  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2"),
    title = "Row Z-Score",
    border = "black",
    title_position = "leftcenter-rot"))
ht_2




group_split <- hgps_correct_order$Individual      

group_split <- factor(group_split, levels = c("NEO", "HGPS 169"))


ht <- Heatmap(
  vsd_zscore_goi,
  cluster_columns = F,
  column_split = group_split,
  top_annotation = ha1,
  col = col_fun,
  border_gp = gpar(col = "black", lty = 1),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 90,
  row_title = "SenMayo DEGs in HGPS vs NEO",
  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2"),
    title = "Row Z-Score",
    border = "black",
    title_position = "leftcenter-rot"))
ht


## SenMayo DEGs in NEO + DS Full Overlap ---------------------------------------

vsd <- vst(dds_neo_ds_umk57, blind=TRUE) ## Blind to ignore conditions

degs_ds_vs_neo <- df_ds_neo %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::pull(gene)

degs_ds <- read.csv("output_data\\excel_results\\deseq2_ds_umk57_vs_dmso.csv") %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::pull(gene)

senamyo_degs_ds_vs_neo <- intersect(
  degs_ds_vs_neo,
  senmayo_genes
)

senamyo_degs_ds_vs_neo <- intersect(
  senamyo_degs_ds_vs_neo,
  degs_ds
)

## Create heatmap with genes of interest:

vsd_df <- assay(vsd)
vsd_zscore <- scale(t(vsd_df)) %>% as.data.frame()

ds_correct_order <- sample_info_correct %>% 
  dplyr::filter(Disease %in% c("Healthy", "DS")) %>% 
  dplyr::arrange(Treatment, desc(Individual))

vsd_zscore <- vsd_zscore[row.names(ds_correct_order),]  

vsd_zscore_goi <- vsd_zscore %>% as.data.frame()
vsd_zscore_goi <- vsd_zscore_goi[,senamyo_degs_ds_vs_neo]


col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-2, 2)) ## Necessary to assign values to colors

ha1 = rowAnnotation(Individual = ds_correct_order$Individual,
                    Treatment = ds_correct_order$Treatment,
                    col = list(
                      Individual = c(
                        "NEO" = "gray",
                        "14y_DS" = "blue",
                        "5y_DS" = "magenta"
                      ),
                      Treatment = c(
                        "DMSO" = "gray",
                        "UMK57" = "cyan"
                      )
                    ))



group_split <- ds_correct_order$Disease      

group_split <- factor(group_split, levels = c("Healthy", "DS"))


ht <- Heatmap(
  vsd_zscore_goi,
  cluster_columns = T,
  cluster_rows = F,
  show_row_names = F,
  row_split = group_split,
  column_km = 2,
  left_annotation = ha1,
  col = col_fun,
  border_gp = gpar(col = "black", lty = 1),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 90,
  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2"),
    title = "Row Z-Score",
    border = "black",
    title_position = "leftcenter-rot"))
ht

ht_2 + ht

## Save Manually 13x3

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
