#' Siwei rewrite 18 Oct 2023
#' SplitNCigarBam
#' Since the bams output from bowtie2 pipeline has been added for readgroups,
#' sorted, and indexed, skip the first three steps
#'
#'
#' Create the VCF (Variant Call Format) file
#' @param Directory The Path of to the BAM Directory
#' @param Genome_Fa Path for whole genome FASTQ file
#' @param Picard_Path Path to the Picard directory, with all the JAR files
#' @param GATK_Path Path to the GATK directory, where the GenomeAnalysisTK.jar file is located
#' @param bam_file_name input bam name
#' @param temp_dir temp directory
#' @export
#' @return None


SplitNCigarBam <-
  function(Directory,
           Genome_Fa,
           Picard_Path,
           GATK_Path,
           bam_file_name,
           temp_dir) {

    print(getwd())
    print(temp_dir)

    # setwd(Directory)

    # create temp_dir if non-existent
    if (!dir.exists(temp_dir)) {
      dir.create(temp_dir)
    } else {
      unlink(temp_dir,
             recursive = T)
      dir.create(temp_dir)
    }
    # comm = paste("java -jar ",
    #              Picard_Path,
    #              "picard.jar AddOrReplaceReadGroups  I=",
    #              Directory,
    #              "accepted_hits.bam
    #              O=accepted_hits_rg.bam ID=\"n\" LB=\"lb\" PL=\"illumina\" PU=\"pu\" SM=\"ES\"",
    #              sep="")
    # print("======================")
    # print(comm)
    # print("======================")
    # system(comm)
    # setwd("~/")
    #
    # comm = paste("java -jar ",Picard_Path, "picard.jar ReorderSam I=",Directory, "accepted_hits_rg.bam O=",Directory, "accepted_hits_rg_sorted.bam R=",Genome_Fa, sep="")
    # print("======================")
    # print(comm)
    # print("======================")
    # system(comm)
    #
    # comm = paste("java -jar ",Picard_Path, "picard.jar BuildBamIndex I=",Directory, "accepted_hits_rg_sorted.bam", sep="")
    # print("======================")
    # print(comm)
    # print("======================")
    # system(comm)

    comm = paste(GATK_Path,
                 '/gatk ',
                 "SplitNCigarReads -R ",
                 Genome_Fa,
                 " -I ",
                 bam_file_name,
                 " -O ",
                 temp_dir,
                 '/',
                 str_replace_all(string = basename(bam_file_name),
                                 pattern = "\\.bam$",
                                 replacement = "\\_split.bam"),
                 sep = "")
    altComm = paste(comm,
                    "--fix_misencoded_quality_scores")
    print("======================")
    print(comm)
    print("======================")
    system(comm)
    # if (!file.exists(paste0(dirname(bam_file_name),
    #                         '/',
    #                         temp_dir,
    #                         str_replace_all(string = basename(bam_file_name),
    #                                         pattern = "\\.bam$",
    #                                         replacement = "\\_split.bam")))) {
    #   print("======================")
    #   print("Alternative Command")
    #   print(altComm)
    #   print("======================")
    #   system(altComm)
    #
    # }


    # comm = paste(GATK_Path,
    #              '/gatk ',
    #              "HaplotypeCaller -R ",
    #              Genome_Fa,
    #              " -I ",
    #              temp_dir,
    #              '/',
    #              str_replace_all(string = basename(bam_file_name),
    #                              pattern = "\\.bam$",
    #                              replacement = "\\_split.bam"),
    #              " --dont-use-soft-clipped-bases true ",
    #              "-stand-call-conf 10.0 ",
    #              "--native-pair-hmm-threads 24 ",
    #              "-O ",
    #              Directory,
    #              '/',
    #              "variants_output.vcf",
    #              sep = "")
    # print("======================")
    # print(comm)
    # print("======================")
    # system(comm)
    #
    # cleanup
    # unlink(x = temp_dir,
    #        recursive = T)

  }
