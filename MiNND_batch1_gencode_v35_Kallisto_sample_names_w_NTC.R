# Siwei 31 May 2024
# load kallisto-quantified per-transcript data

# init ####
{
  library(readr)
  library(readxl)
  library(edgeR)
  library(stringr)

  library(tximport)

  library(GenomicFeatures)
}


# GenomicFeatures::makeTxDbFromEnsembl(organism = "Homo sapiens",
#                                      release = 108)
txdb_Ensembl_108 <-
  GenomicFeatures::makeTxDbFromGFF(file = "/home/zhangs3/Data/Databases/Genomes/hg38/Homo_sapiens.chr.GRCh38.108.chr.gtf",
                                   organism = "Homa sapiens",
                                   taxonomyId = 9606)
columns(txdb_Ensembl_108)
keys(txdb_Ensembl_108)

head(keys(txdb_Ensembl_108,
          keytype = "TXNAME"))

load("~/backuped_space/Siwei_misc_R_projects/Alena_RNASeq_23Aug2023/ENSG_gene_index.RData")

ENSG_anno_gene_indexed[ENSG_anno_gene_indexed$Gene_Symbol %in% c("TRIO",
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
                                                                 "SCN2A"), 1]

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

files_import <-
  dir(path = "/home/zhangs3/Data/FASTQ/RNASeq_Novogene_10Oct2023/fastq_kallisto/old_kallisto_output",
      pattern = "*.tsv",
      full.names = T,
      recursive = T)

sample_names <-
  str_split(string = files_import,
            pattern = "\\/",
            simplify = T)
sample_names <-
  unlist(sample_names[, 9])

# 07 Jun: use NTC
sample_sans_NTC <- sample_names
#   sample_names[!(sample_names %in% "V3865")]
#
# files_import <- files_import
#   files_import[!(sample_names %in% "V3865")]

tx.kallisto <-
  tximport(files_import,
           type = "kallisto",
           txIn = T,
           txOut = T,
           countsFromAbundance = "scaledTPM",
           txIdCol = "target_id")

sum(str_detect(rownames(tx.kallisto$counts),
               pattern = "ENST00000617797"))

names(files_import) <- sample_sans_NTC
# save.image("MiNND_batch1_kallisto_tx.RData")
# load("MiNND_batch1_kallisto_tx.RData")

ENSG_genes_list <-
  ENSG_anno_gene_indexed[ENSG_anno_gene_indexed$Gene_Symbol %in% KD_genes_list, 1]
ENSG_genes_list$Gene_symbol <-
  ENSG_anno_gene_indexed[ENSG_anno_gene_indexed$Gene_Symbol %in% KD_genes_list, 2]
ENSG_genes_list$Gene_symbol <-
  unlist(ENSG_genes_list$Gene_symbol)

sample_gene_lookup <-
  data.frame(Gene_symbol = KD_genes_list,
             Sample_id = sample_sans_NTC)
# sample_gene_lookup$Sample_id <-
#   unlist(sample_gene_lookup$Sample_id)

ENSG_gene_sample_lookup <-
  base::merge(x = ENSG_genes_list,
              y = sample_gene_lookup,
              by.x = "Gene_symbol",
              by.y = "Gene_symbol")
ENSG_gene_sample_lookup <-
  ENSG_gene_sample_lookup[order(ENSG_gene_sample_lookup$Sample_id), ]

names(tx.kallisto$abundance)

tx.abundance <-
  tx.kallisto$abundance
nrow(tx.abundance)
rownames(tx.abundance) <-
  str_split(string = rownames(tx.abundance),
            pattern = "\\.",
            simplify = T)[, 1]

sum(str_detect(rownames(tx.abundance),
               pattern = "ENST00000617797"))

# tx.abundance <-
#   as.data.frame(t(tx.abundance))

gene2transcripts.list <-
  vector(mode = "list",
         length = length(KD_genes_list))
names(gene2transcripts.list) <-
  ENSG_gene_sample_lookup$Geneid

tx.extracted.kallisto <-
  vector(mode = "list",
         length = length(KD_genes_list))
names(tx.extracted.kallisto)

for (i in 1:length(gene2transcripts.list)) {
  print(i)
  # names(gene2transcripts.list)[i] <-
  #   ENSG_gene_sample_lookup$Gene_symbol[i]
  returned_transcript_vector <-
    unlist(mapIds(x = txdb_Ensembl_108,
                  keys = ENSG_gene_sample_lookup$Geneid[i],
                  keytype = "GENEID",
                  column = "TXNAME",
                  multiVals = "list"))
  names(tx.extracted.kallisto)[i] <-
    ENSG_gene_sample_lookup$Gene_symbol[i]
  tx.extracted.kallisto.value <-
    as.data.frame(tx.abundance[rownames(tx.abundance) %in% returned_transcript_vector, ])
  if (ncol(tx.extracted.kallisto.value) == 1) {
    tx.extracted.kallisto.value <-
      as.data.frame(t(tx.extracted.kallisto.value))
  }

  tx.extracted.kallisto[[i]] <-
    tx.extracted.kallisto.value
  # gene2transcripts.list[[i]] <-
  #   unlist(mapIds(x = txdb_Ensembl_108,
  #                 keys = ENSG_gene_sample_lookup$Geneid[i],
  #                 keytype = "GENEID",
  #                 column = "TXNAME",
  #                 multiVals = "list"))
  # print(all_transcripts_list)
}

nrow(tx.extracted.kallisto[[4]])
ncol(tx.extracted.kallisto[[4]])


for (i in 1:length(tx.extracted.kallisto)) {
  print(i)
  assembled_rownames <-
    str_c(rownames(tx.extracted.kallisto[[i]]),
          names(tx.extracted.kallisto)[i],
          sep = "-")
  print(paste("rownames length =",
              length(assembled_rownames)))
  assembled_colnames <-
    str_c(ENSG_gene_sample_lookup$Sample_id,
          ENSG_gene_sample_lookup$Gene_symbol,
          sep = "-")
  print(paste("colnames length =",
              length(assembled_colnames)))
  rownames(tx.extracted.kallisto[[i]]) <-
    assembled_rownames
  colnames(tx.extracted.kallisto[[i]]) <-
    assembled_colnames
}

master.tx.table.output <-
  do.call(what = rbind,
          args = tx.extracted.kallisto)
rownames(master.tx.table.output) <-
  str_split(string = rownames(master.tx.table.output),
            pattern = "-",
            simplify = T)[, 1]

write.table(master.tx.table.output,
            file = "master_tx_table.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

mapIds(x = txdb_Ensembl_108,
       keys = "ENSG00000023516",
       keytype = "GENEID",
       column = "TXNAME")


transcript_list_kallisto <-
  readxl::read_excel(path = "master_tx_table_gencode_v35_rev.xlsx",
                     sheet = 1)
# transcript_ENSG <-
#   str_split(string = rownames(transcript_list_kallisto),
#             pattern = "\\.ENST",
#             simplify = T)[, 2]


tx_selected_transcripts_kallisto <-
  # tx.kallisto$
  # tx.abundance <-
  tx.kallisto$abundance
nrow(tx_selected_transcripts_kallisto)
nrow(transcript_list_kallisto)

rownames(tx_selected_transcripts_kallisto)
tx_selected_transcripts_kallisto_rownames <-
  rownames(tx_selected_transcripts_kallisto)
tx_selected_transcripts_kallisto_rownames <-
  str_split(rownames(tx_selected_transcripts_kallisto),
            pattern = "\\|",
            simplify = T)[, 1]
colnames(tx_selected_transcripts_kallisto) <-
  sample_names

sum(duplicated(rownames(transcript_list_kallisto)))
sum(duplicated(str_split(rownames(tx_selected_transcripts_kallisto),
                         pattern = "\\.",
                         simplify = T)[, 1]))
sum(duplicated(str_split(transcript_list_kallisto$Gene.TxID,
                         pattern = "\\.",
                         simplify = T)[, 2]))

dup_detect <-
  transcript_list_kallisto[duplicated(str_split(rownames(transcript_list_kallisto),
                                                pattern = "\\.",
                                                simplify = T)[, 2]), ]

sum(str_detect(string = rownames(tx_selected_transcripts_kallisto),
               pattern = "ENST00000261917"))
sum(str_detect(string = transcript_list_kallisto$Gene.TxID,
               pattern = "ENST00000261917"))

sum(str_split(rownames(tx_selected_transcripts_kallisto),
              pattern = "\\.",
              simplify = T)[, 1] %in%
      str_split(transcript_list_kallisto$Gene.TxID,
                pattern = "\\.",
                simplify = T)[, 2])

imported_tx_list <-
  data.frame(TxID = str_split(rownames(tx_selected_transcripts_kallisto),
                              pattern = "\\.",
                              simplify = T)[, 1])
imported_tx_list$dup_mark <-
  "non-dup"
xlsx_tx_list <-
  data.frame(TxID = str_split(transcript_list_kallisto$Gene.TxID,
                              pattern = "\\.",
                              simplify = T)[, 2])

merged_dup_list <-
  base::merge(x = xlsx_tx_list,
              y = imported_tx_list,
              by.x = "TxID",
              by.y = "TxID",
              all.x = T)

rownames(tx_selected_transcripts_kallisto)[str_detect(string = rownames(tx_selected_transcripts_kallisto),
                                                      pattern = "ENST00000617797")]


tx_selected_transcripts_kallisto_output <-
  tx_selected_transcripts_kallisto[str_split(rownames(tx_selected_transcripts_kallisto),
                                             pattern = "\\.",
                                             simplify = T)[, 1] %in%
                                     str_split(transcript_list_kallisto$Gene.TxID,
                                               pattern = "\\.",
                                               simplify = T)[, 2], ]
colnames(tx_selected_transcripts_kallisto)

nrow(tx_selected_transcripts_kallisto)
nrow(tx_selected_transcripts_kallisto_output)
length(tx_selected_transcripts_kallisto[, colnames(tx_selected_transcripts_kallisto) %in% "V3865"])


NTC_transcripts_kallisto <-
  data.frame(V3865.NTC = tx_selected_transcripts_kallisto_output[, colnames(tx_selected_transcripts_kallisto) %in% "V3865"])
rownames(NTC_transcripts_kallisto) <-
  str_split(rownames(NTC_transcripts_kallisto),
            pattern = "\\|",
            simplify = T)[, 1]
rownames(NTC_transcripts_kallisto) <-
  str_split(rownames(NTC_transcripts_kallisto),
            pattern = "\\.",
            simplify = T)[, 1]
NTC_transcripts_kallisto$TxID <-
  rownames(NTC_transcripts_kallisto)

transcript_list_kallisto$TxID <-
  str_split(transcript_list_kallisto$Gene.TxID,
            pattern = "\\.",
            simplify = T)[, 2]

transcript_all_w_NTC <-
  base::merge(x = transcript_list_kallisto,
              y = NTC_transcripts_kallisto,
              by = "TxID")
transcript_all_w_NTC <-
  base::merge(x = NTC_transcripts_kallisto,
              y = transcript_list_kallisto,
              by = "TxID")

transcript_all_w_NTC <-
  transcript_all_w_NTC[transcript_all_w_NTC$TxID %in%
                         transcript_list_kallisto$TxID, ]
transcript_all_w_NTC <-
  transcript_all_w_NTC[transcript_list_kallisto$TxID %in%
                         transcript_all_w_NTC$TxID, ]

transcript_all_w_NTC <-
  transcript_all_w_NTC[base::match(x = transcript_list_kallisto$TxID,
                                   table = transcript_all_w_NTC$TxID), ]

# write.table()

write.table(transcript_all_w_NTC,
            file = "master_tx_table_gencode_v35_w_NTC.tsv",
            quote = F, sep = "\t",
            row.names = F, col.names = T)
