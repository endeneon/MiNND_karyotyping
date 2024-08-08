# 19 Oct 2023 Siwei
# Write a loop to walk through all ready-made BAM files
# SplitCigarNReads skipped, should not cause major problem

# 06 Jan 2024 Siwei
# add a few lines so the code can run independently (i.e. multiple instances
# can run at the same time) in different directories, in parallel)

# Libraries
{
  library(zoo)
  library(gplots)

  library(stringr)
}

setwd(".")
print(paste("working directory =",
            getwd()))

# ! verify if the directory storing input BAMs exists ! ####
bam_file_dir = "." # ! manully set bam path here !!!
if (!dir.exists(bam_file_dir)) {
  print("BAM input dir is not set!")
  q(save = "no",
    status = 1)
}
print(paste("bam_file_dir =",
            bam_file_dir))

# init environment ####
Picard_Path <- "/home/zhangs3/Data/Tools/gatk-4.1.8.1"
Genome_Fa <- "/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
GATK_Path <- "/home/zhangs3/Data/Tools/gatk-4.1.8.1"
Organism <- "Human"

temp_dir <- "./temp_dir"
# test if temp_dir exists
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir)
}

# source all functions ####
function_list_to_source <-
  list.files(path = "/backup/Siwei_misc_R_projects/eSNPKaryotyping/eSNPKaryotyping/R",
             pattern = "*.R",
             full.names = T)
for (i in 1:length(function_list_to_source)) {
  source(function_list_to_source[i])
}
# make a eSNPKaryotyping-compatible dbSNP file (one-off) #####
# ! Use $Databases/Genomes/hg38/common_snp150_by_chr ! #####
# Edit_dbSNP_Files(Directory = "/home/zhangs3/Data/Databases/Genomes/hg38/common_snp150_by_chr",
#                  File_Name = "Edited_common_chr",
#                  Organism = "Human")

# dir.create("working")
# dir.create("plot_output")


bam_file_list <-
  list.files(path = bam_file_dir,
             full.names = T,
             recursive = F,
             pattern = "*bam$")
print(bam_file_list)

Directory <-
  "./working"
if (!dir.exists(Directory)) {
  dir.create(Directory)
}

plot_output_directory <-
  "./plot_output"
if (!dir.exists(plot_output_directory)) {
  dir.create(plot_output_directory)
}


try({
  dev.off()
})

k <- 1
for (k in 1:length(bam_file_list)) {
  print(paste(k, bam_file_list[k], sep = " "))
  print(date())
  # set the file name and title for the PDF
  pdf_series_name <-
    basename(bam_file_list[k])
  pdf_series_name <-
    gsub(pattern = ".bam",
         replacement = "",
         x = pdf_series_name)
  print(pdf_series_name)

  # need to pass temp_dir every time
  CreateVCF(Directory = Directory,
            Genome_Fa = Genome_Fa,
            Picard_Path = Picard_Path,
            GATK_Path = GATK_Path,
            bam_file_name = bam_file_list[k],
            temp_dir = temp_dir)

  EditVCF(Directory = Directory,
          Organism = Organism,
          temp_dir = temp_dir)

  # read the table in
  VCF_table <-
    read.delim(file = paste(Directory,
                            "variantTable.csv",
                            sep = "/"),
               header = T, sep = "\t",
               quote = "",
               dec = ".")
  VCF_table$chr <-
    as.numeric(VCF_table$chr)
  VCF_table <-
    VCF_table[order(VCF_table$chr,
                    VCF_table$position), ]
  VCF_table <-
    VCF_table[VCF_table$chr > 0, ]
  #### backup
  # VCF_table_backup <- VCF_table
  # VCF_table <- VCF_table_backup
  # VCF_table$chr <- paste("chr", VCF_table$chr, sep = "")


  ###### Return MajorMinorCalc results
  MajorMinorCal_results <-
    MajorMinorCalc(Table = VCF_table,
                   minDP = 20, # ! change back to 20 later !
                   maxDP = 1000,
                   minAF = 0.2)

  ##### Plot Allelic ratio along the genome for duplication detection
  # setwd("../")
  print(getwd())
  if (!dir.exists(plot_output_directory)) {
    dir.create(plot_output_directory)
  }
  pdf(file = paste(plot_output_directory,
                   '/',
                   pdf_series_name,
                   "_genome.pdf",
                   sep = ""),
      paper = "USr",
      title = pdf_series_name)

  PlotGenome(orderedTable = MajorMinorCal_results,
             Window = 151,
             Ylim = 3,
             PValue = TRUE,
             Organism = Organism)
  try({
    dev.off()
  })
  # save current working directory
  # current_working_directory <- getwd()
  #pdf(plot)
  # setwd(plot_output_directory) # send plot to plot output directory

  # setwd(current_working_directory)

  ##### intersect the observed SNP with common SNP from dbSNPs
  # Siwei: the Genome_Fa_dict is simply a waste of memory...
  # moreover, it did not consider the situation that if the .dict has
  # additional header but not started with @SQ
  # setwd("~/backuped_space/Siwei_misc_R_projects/eSNPKaryotyping/eSNPKaryotyping")
  # setwd(".")
  # print(getwd())

  tbl_DeletionTable_output <-
    DeletionTable(Directory = Directory,
                  Table = MajorMinorCal_results,
                  dbSNP_Data_Directory = "/home/zhangs3/Data/Databases/Genomes/hg38/common_snp150_by_chr",
                  dbSNP_File_Name = "Edited_Common_chr",
                  Genome_Fa_dict = Genome_Fa,
                  Organism = "Human",
                  source_bam = bam_file_list[k],
                  temp_dir = temp_dir)


  # 9. Plot each SNP, without any summarization
  # ! very time-consuming !18Oct2
  # setwd("../")
  # pdf(file = paste(plot_output_directory,
  #                  '/',
  #                  pdf_series_name,
  #                  "_zygosity_single.pdf",
  #                  sep = ""),
  #     paper = "USr")
  # Plot_Zygosity_Sinle(Table = tbl_DeletionTable_output,
  #                     Organism = "Human")
  # # Argument: 1. Table - The LOH table containing the output of the DeletionTable function
  # #           2. Organism - "Human" or "Mouse"
  # dev.off()




  # 10. Plot blocks of heterozygous and homozygous SNPs
  print(getwd())
  pdf(file = paste(plot_output_directory,
                   '/',
                   pdf_series_name,
                   "_zygosity_blocks.pdf", sep = ""),
      paper = "USr",
      title = pdf_series_name)
  Plot_Zygosity_Blocks(Table = tbl_DeletionTable_output,
                       Window = 1500000,
                       Max = 6,
                       Max2 = 20,
                       Organism = "Human")
  # Argument: 1. Table - The deletion table containing the output of the DeletionTable function
  #           2. window - the block size in bp, usually 1500000
  #           3. Max - How many Heterozygouse SNP need to be in a block to get the full color, usually 6
  #           4. Max2 - How many Homozygouse SNP need to be in a block to get the full color, usually 60
  #           5. Organism - "Human" or "Mouse"
  try({
    dev.off()
  })
}

###################
# bam_file_list[1]
# paste("I=",Directory,bam_file_list[1])

# Edit_dbSNP_Files(Directory = "~/1TB/Databases/hg38/common_snp150_by_chr/",
#                  File_Name = "dbSNP_150_common_chr",
#                  Organism = "Human")

## read hg38 chr sizes
# hg38_chr_size_info <-
#   read.table(file = "/home/zhangs3/Data/Databases/Genomes/hg38/GRCh38p13.primary_assembly.genome.dict",
#              header = F,
#              sep = "\t",
#              skip = 1)
# hg38_chr_size_info$length <-
#   str_split(string = hg38_chr_size_info$V3,
#             pattern = ":",
#             simplify = T)[, 2]
# # as.numeric(hg38_chr_size_info$length[1:24])
#
# ## read hg38 centromere sizes
# hg38_chr_cen_info <-
#   read.table(file = "/home/zhangs3/Data/Databases/Genomes/hg38/hg38_centromere_positions.bed",
#              header = F,
#              sep = "\t",
#              skip = 0)
# hg38_chr_cen_info <-
#   hg38_chr_cen_info[(1:nrow(hg38_chr_cen_info)) %% 2 == 0, ]
# # as.numeric(hg38_chr_cen_info$V2)
#
#
# # need the split NCigar Bam...
# SplitNCigarBam(Directory = Directory,
#                Genome_Fa = Genome_Fa,
#                Picard_Path = Picard_Path,
#                GATK_Path = GATK_Path,
#                bam_file_name = bam_file_list[k],
#                temp_dir = temp_dir)
