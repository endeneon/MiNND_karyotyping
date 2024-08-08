# Siwei 07Jun 2024
# plot fig. 5a

# init ####
{
  library(gplots)
  library(RColorBrewer)
  library(readxl)

  library(readr)
  library(edgeR)

  library(stringr)
}

load("~/backuped_space/Siwei_misc_R_projects/Alena_RNASeq_23Aug2023/ENSG_gene_index.RData")
KD_genes_list <-
  rev(c("TRIO",
        "XPO7",
        "KDM6B",
        "ANKRD11",
        "SETD1A",
        "CUL1",
        "AKAP11",
        "GRIA3",
        "GRIN2A",
        "KMT2C",
        "GABRA1",
        "ARID1B",
        "DLL1",
        "SP4",
        "RB1CC1",
        "SHANK3",
        "ASH1L",
        "SMARCC2",
        "HCN4",
        "CACNA1G",
        "CHD8",
        "SCN2A"))

b1_metadata_with_NTC <-
  read_csv("~/Data/Alexi_Duhe/miind_scratchpaper/all_batches_read_tables/b1_metadata_with_NTC.csv")

df_raw <-
  read.table(file = "~/Data/Alexi_Duhe/miind_scratchpaper/all_batches_read_tables/batch_1_raw_reads.txt",
             header = T,
             sep = "\t")

df_2_DGE <-
  df_raw
df_2_DGE$Geneid <-
  str_split(df_2_DGE$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
colnames(df_2_DGE) <-
  str_remove_all(colnames(df_2_DGE),
                 pattern = "_")

df_2_DGE <-
  df_2_DGE[!duplicated(df_2_DGE$Geneid), ]
df_2_DGE <-
  merge(x = df_2_DGE,
        y = ENSG_anno_gene_indexed,
        by.x = "Geneid",
        by.y = "Geneid")

df_2_DGE <-
  df_2_DGE[!duplicated(df_2_DGE$Gene_Symbol), ]

df_gene_list <-
  df_2_DGE$GeneSymbol

df_2_DGE$Geneid <- NULL
df_2_DGE$GeneSymbol <- NULL

rownames(df_2_DGE) <-
  df_gene_list

df_2_DGE <-
  DGEList(counts = as.matrix(df_2_DGE),
          samples = colnames(df_2_DGE),
          genes = rownames(df_2_DGE),
          remove.zeros = T)

df_cpm <-
  cpm(df_2_DGE,
      normalized.lib.sizes = T,
      log = F)

df_trimmed_cpm <-
  as.data.frame(df_cpm[rownames(df_cpm) %in% KD_genes_list, ])

df_cpm_from_xlsx <-
  read_excel(path = "../FACS/Fig_5a.xlsx")
