#' @title ascat.loadData
#' @description Function to read in SNP array data
#' @details germline data files can be NULL - in that case these are not read in
#' @param Tumor_LogR_file file containing logR of tumour sample(s)
#' @param Tumor_BAF_file file containing BAF of tumour sample(s)
#' @param Germline_LogR_file file containing logR of germline sample(s), NULL
#' @param Germline_BAF_file file containing BAF of germline sample(s), NULL
#' @param chrs a vector containing the names for the chromosomes (e.g. c(1:22, "X"))
#' @param gender a vector of gender for each cases ("XX" or "XY"). Default = all female ("XX")
#' @param sexchromosomes a vector containing the names for the sex chromosomes. Default = c("X", "Y")
#' @param genomeVersion a string (either 'hg19' or 'hg38') so nonPAR coordinates on X can be stored, NULL
#' @param isTargetedSeq a boolean indicating whether data come from a targeted sequencing experiment. Default = F
#'
#' @return ascat data structure containing:\cr
#' 1. Tumor_LogR data matrix\cr
#' 2. Tumor_BAF data matrix\cr
#' 3. Tumor_LogR_segmented: placeholder, NULL\cr
#' 4. Tumor_BAF_segmented: placeholder, NULL\cr
#' 5. Germline_LogR data matrix\cr
#' 6. Germline_BAF data matrix\cr
#' 7. SNPpos: position of all SNPs\cr
#' 8. ch: a list containing vectors with the indices for each chromosome (e.g. Tumor_LogR[ch[[13]], ] will output the Tumor_LogR data of chromosome 13\cr
#' 9. chr: a list containing vectors with the indices for each distinct part that can be segmented separately (e.g. chromosome arm, stretch of DNA between gaps in the array design)\cr
#' 10. chrs: a vector containing chromosome names\cr
#' 11. samples: a vector containing sample name(s)\cr
#' 12. gender: a vector of gender for each cases ("XX" or "XY"). Default = NULL: all female ("XX")\cr
#' 13. sexchromosomes: a vector containingg names of sex chromosomes\cr
#' 14. X_nonPAR: a vector of two values (start and stop) to define where the nonPAR region is on X\cr
#' 15. isTargetedSeq:  boolean indicating whether data come from a targeted sequencing experiment\cr
#' 16. failedarrays: placeholder, NULL\cr
#'
#' @export
#'
ascat.loadData = function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file = NULL, Germline_BAF_file = NULL, chrs = c(1:22, "X", "Y"), gender = NULL, sexchromosomes = c("X", "Y"), genomeVersion=NULL, isTargetedSeq=FALSE) {

  stopifnot(length(isTargetedSeq)==1 && isTargetedSeq %in% c(TRUE, FALSE))
  # read in SNP array data files
  print.noquote("Reading Tumor LogR data...")
  Tumor_LogR <- read.table(Tumor_LogR_file, header=TRUE, row.names=1, comment.char="", sep = "\t", check.names=FALSE)
  print.noquote("Reading Tumor BAF data...")
  Tumor_BAF <- read.table(Tumor_BAF_file, header=TRUE, row.names=1, comment.char="", sep = "\t", check.names=FALSE)

  #infinite values are a problem - change those
  Tumor_LogR[Tumor_LogR==-Inf]=NA
  Tumor_LogR[Tumor_LogR==Inf]=NA

  Germline_LogR = NULL
  Germline_BAF = NULL
  if (!is.null(Germline_LogR_file)) {
    print.noquote("Reading Germline LogR data...")
    Germline_LogR <- read.table(Germline_LogR_file, header=TRUE, row.names=1, comment.char="", sep = "\t", check.names=FALSE)
    print.noquote("Reading Germline BAF data...")
    Germline_BAF <- read.table(Germline_BAF_file, header=TRUE, row.names=1, comment.char="", sep = "\t", check.names=FALSE)

    #infinite values are a problem - change those
    Germline_LogR[Germline_LogR==-Inf]=NA
    Germline_LogR[Germline_LogR==Inf]=NA
  }

  # make SNPpos vector that contains genomic position for all SNPs and remove all data not on chromosome 1-22, X, Y (or whatever is given in the input value of chrs)
  print.noquote("Registering SNP locations...")
  SNPpos <- Tumor_LogR[, 1:2]
  SNPpos = SNPpos[SNPpos[, 1] %in% chrs, ]

  # if some chromosomes have no data, just remove them
  chrs = intersect(chrs, unique(SNPpos[, 1]))

  Tumor_LogR = Tumor_LogR[rownames(SNPpos), c(-1, -2), drop=FALSE]
  Tumor_BAF = Tumor_BAF[rownames(SNPpos), c(-1, -2), drop=FALSE]
  # make sure it is all converted to numerical values
  for (cc in 1:dim(Tumor_LogR)[2]) {
    Tumor_LogR[, cc]=as.numeric(as.vector(Tumor_LogR[, cc]))
    Tumor_BAF[, cc]=as.numeric(as.vector(Tumor_BAF[, cc]))
  }

  if (!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[rownames(SNPpos), c(-1, -2), drop=FALSE]
    Germline_BAF = Germline_BAF[rownames(SNPpos), c(-1, -2), drop=FALSE]
    for (cc in 1:dim(Germline_LogR)[2]) {
      Germline_LogR[, cc]=as.numeric(as.vector(Germline_LogR[, cc]))
      Germline_BAF[, cc]=as.numeric(as.vector(Germline_BAF[, cc]))
    }
  }

  # sort all data by genomic position
  last = 0
  ch = list()
  SNPorder = vector(length=dim(SNPpos)[1])
  for (i in 1:length(chrs)) {
    chrke = SNPpos[SNPpos[, 1]==chrs[i], ]
    chrpos = chrke[, 2]
    names(chrpos) = rownames(chrke)
    chrpos = sort(chrpos)
    ch[[i]] = (last+1):(last+length(chrpos))
    SNPorder[ch[[i]]] = names(chrpos)
    last = last+length(chrpos)
  }
  SNPpos = SNPpos[SNPorder, ]
  Tumor_LogR=Tumor_LogR[SNPorder, , drop=FALSE]
  Tumor_BAF=Tumor_BAF[SNPorder, , drop=FALSE]

  if (!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[SNPorder, , drop=FALSE]
    Germline_BAF = Germline_BAF[SNPorder, , drop=FALSE]
  }

  # split the genome into distinct parts to be used for segmentation (e.g. chromosome arms, parts of genome between gaps in array design)
  if (!isTargetedSeq) {
    print.noquote("Splitting genome in distinct chunks...")
    chr = split_genome(SNPpos)
  } else {
    chr=ch
  }

  if (is.null(gender)) {
    gender = rep("XX", dim(Tumor_LogR)[2])
  }

  if (!is.null(genomeVersion)) {
    if (genomeVersion=="hg19") {
      X_nonPAR=c(2699521, 154931043)
    } else if (genomeVersion=="hg38") {
      X_nonPAR=c(2781480, 155701382)
    } else {
      stop("genomeVersion must be either \'hg19\' or \'hg38\'.")
    }
  } else {
    X_nonPAR=NULL
  }

  return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF,
              Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL,
              Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF,
              SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs,
              samples = colnames(Tumor_LogR), gender = gender,
              sexchromosomes = sexchromosomes, X_nonPAR = X_nonPAR,
              isTargetedSeq = isTargetedSeq, failedarrays = NULL))
}

# helper function to split the genome into parts
split_genome = function(SNPpos) {

  # look for gaps of more than 5Mb (arbitrary treshold to account for big centremeres or other gaps) and chromosome borders
  bigHoles = which(diff(SNPpos[, 2])>=5000000)+1
  chrBorders = which(SNPpos[1:(dim(SNPpos)[1]-1), 1]!=SNPpos[2:(dim(SNPpos)[1]), 1])+1

  holes = unique(sort(c(bigHoles, chrBorders)))

  # find which segments are too small
  #joincandidates=which(diff(c(0, holes, dim(SNPpos)[1]))<200)

  # if it's the first or last segment, just join to the one next to it, irrespective of chromosome and positions
  #while (1 %in% joincandidates) {
  #  holes=holes[-1]
  #  joincandidates=which(diff(c(0, holes, dim(SNPpos)[1]))<200)
  #}
  #while ((length(holes)+1) %in% joincandidates) {
  #  holes=holes[-length(holes)]
  #  joincandidates=which(diff(c(0, holes, dim(SNPpos)[1]))<200)
  #}

  #while(length(joincandidates)!=0) {
  # the while loop is because after joining, segments may still be too small..

  #startseg = c(1, holes)
  #endseg = c(holes-1, dim(SNPpos)[1])

  # for each segment that is too short, see if it has the same chromosome as the segments before and after
  # the next always works because neither the first or the last segment is in joincandidates now
  #previoussamechr = SNPpos[endseg[joincandidates-1], 1]==SNPpos[startseg[joincandidates], 1]
  #nextsamechr = SNPpos[endseg[joincandidates], 1]==SNPpos[startseg[joincandidates+1], 1]

  #distanceprevious = SNPpos[startseg[joincandidates], 2]-SNPpos[endseg[joincandidates-1], 2]
  #distancenext = SNPpos[startseg[joincandidates+1], 2]-SNPpos[endseg[joincandidates], 2]

  # if both the same, decide based on distance, otherwise if one the same, take the other, if none, just take one.
  #joins = ifelse(previoussamechr&nextsamechr,
  #               ifelse(distanceprevious>distancenext, joincandidates, joincandidates-1),
  #               ifelse(nextsamechr, joincandidates, joincandidates-1))

  #holes=holes[-joins]

  #joincandidates=which(diff(c(0, holes, dim(SNPpos)[1]))<200)
  #}
  # if two neighboring segments are selected, this may make bigger segments then absolutely necessary, but I'm sure this is no problem.

  startseg = c(1, holes)
  endseg = c(holes-1, dim(SNPpos)[1])

  chr=list()
  for (i in 1:length(startseg)) {
    chr[[i]]=startseg[i]:endseg[i]
  }

  return(chr)
}