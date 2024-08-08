#' EditVCF
#'
#' Edit the VCF (Variant Call Format) file, Creates file with SNPs data at the BAM directory called variantTable.csv
#' @param Directory The Path of to the BAM Directory, also the VCF file
#' @param Organism "Human" or "Mouse"
#' @param temp_dir "temp directory
#' @export
#' @return None

EditVCF <-
  function(Directory,
           Organism,
           temp_dir) {
  print("Editing VCF File")
  Dir = Directory
  # Dir = "."
  file = "variants_output.vcf"
  path = paste(Dir,
               file,
               sep = "/")
  readData = read.delim(path,
                        as.is = T)
  readData <-
    as.character(readData[-c(1:which(readData == "#CHROM") - 1), 1])
  print(readData[1:30])

  jump = 10
  startChr = 1 + jump
  startPos = 2 + jump
  startInfo = 10 + jump

  len = length(readData)
  print(paste("the length of readData is",
              len))
  # here coverted "chr[num] to [num]"
  chrRegex = "^chr(\\w+)$"
  infoRegex = "^([01])\\/([01]):(\\d+)\\,(\\d+):(\\d+):\\d+:\\d+\\,\\d+\\,\\d+$"

  chrVector =  readData[startChr]
  posVector =  readData[startPos]
  infoVector = readData[startInfo]

  cat('===1===')

  while (startInfo + jump < len) {
    startChr = startChr + jump
    startPos = startPos + jump
    startInfo = startInfo + jump
    chrVector = append(chrVector, readData[startChr])
    posVector = append(posVector, readData[startPos])
    infoVector = append(infoVector, readData[startInfo])
  }

  cat('===2===')
  chrNum <-
    gsub(chrRegex,
         "\\1",
         chrVector)
  if (Organism == "Human") {
    chrNum[chrNum == "X"] = "23"
    chrNum[chrNum == "Y"] = "24"
    }

  if (Organism == "Mouse") {
    chrNum[chrNum == "X"] = "20"
    chrNum[chrNum == "Y"] = "21"
    }

  chrNum <- as.numeric(chrNum)
  Karyotape <-
    10 * abs(as.numeric(gsub(infoRegex,
                             "\\1",
                             infoVector)) -
               as.numeric(gsub(infoRegex,
                               "\\2",
                               infoVector)))
  AD1 <- as.numeric(gsub(infoRegex,
                        "\\3",
                        infoVector))
  AD2 <- as.numeric(gsub(infoRegex,
                        "\\4",
                        infoVector))
  DP <- as.numeric(gsub(infoRegex,
                       "\\5",
                       infoVector))

  posVector <- as.numeric(posVector)

  cat('===3===')

  table <-
    data.frame("chr" = chrNum,
               "position" = posVector,
               "AD1" = AD1,
               "AD2" = AD2,
               "DP" = DP,
               "Karyotape" = Karyotape)
  # strip all columns without proper chr names
  table <-
    table[complete.cases(table[, 1:6]), ]

  cat('===4===')

  fileName <- "variantTable.csv"
  pathToSave <-
    paste(Dir,
          fileName,
          sep = "/")
  # table[is.na(table)] <- 0
  write.table(table,
              file = pathToSave,
              sep = "\t",
              row.names = F, col.names = T,
              quote = F)
  return(0)
}
