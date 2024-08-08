# Siwei 30 Oct 2023
# Count knockdown efficiency of MiNND batch 1 data

# init ####
library(readr)
library(readxl)
library(edgeR)
library(stringr)

load("~/backuped_space/Siwei_misc_R_projects/Alena_RNASeq_23Aug2023/ENSG_gene_index.RData")

# load sample table
sample_table <-
  read_excel("input_table/R_input_table_batch1.xlsx")

# load count matrix
df_raw <-
  read_delim("count_matrix/ReadsPerGene_STAR_batch1.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

df_counts <-
  df_raw
df_counts$Geneid <-
  str_split(string = df_counts$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
df_counts <-
  df_counts[!duplicated(df_counts$Geneid), ]
df_counts <-
  merge(df_counts,
        ENSG_anno_gene_indexed,
        by = "Geneid")
df_counts <-
  df_counts[!duplicated(df_counts$Geneid), ]
df_counts <-
  df_counts[!duplicated(df_counts$Gene_Symbol), ]
gene_list <-
  df_counts$Gene_Symbol
df_counts$Geneid <- NULL
df_counts$Gene_Symbol <- NULL
colnames(df_counts) <-
  str_remove_all(string = colnames(df_counts),
                 pattern = "_")
rownames(df_counts) <-
  gene_list

df_DGE <-
  DGEList(counts = as.matrix(df_counts),
          samples = colnames(df_counts),
          genes = gene_list,
          remove.zeros = T)
sum(rownames(df_DGE) %in% sample_table$Gene_symbol)

df_DGE <-
  calcNormFactors(df_DGE)
cpm_DGE <-
  cpm(df_DGE,
      normalized.lib.sizes = T,
      log = F)
cpm_DGE <-
  as.data.frame(cpm_DGE)
df_sample_index <-
  sample_table[sample_table$Gene_symbol %in% rownames(cpm_DGE), ]
df_sample_index <-
  df_sample_index[!(df_sample_index$Gene_symbol %in% "NTC"), ]
df_sample_index <-
  df_sample_index[, c("Gene_symbol",
                      "barcode")]
df_sample_index$barcode <-
  str_replace(string = df_sample_index$barcode,
              pattern = "^V00087",
              replacement = "V")

cpm_DGE$Gene_symbol <-
  rownames(cpm_DGE)

df_cpm_count_samples <-
  merge(x = df_sample_index,
        y = cpm_DGE,
        by = "Gene_symbol")
df_cpm_count_samples <-
  df_cpm_count_samples[order(df_cpm_count_samples$barcode), ]

write.table(df_cpm_count_samples,
            file = "cpm_count_tables_MiNND_batch1_30Oct2023.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")
