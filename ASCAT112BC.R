# read data: LRR and BAF, both from tumor and blood material;
Tumor_LogR <- read.table("Tumor_LRR.txt", header=T, row.names=1)
Blood_LogR <- read.table("Blood_LRR.txt", header=T, row.names=1)
Tumor_BAllele <- read.table("Tumor_BAF.txt", header=T, row.names=1)
Blood_BAllele <- read.table("Blood_BAF.txt", header=T, row.names=1)

# SNP annotation is in the first two columns (right after the SNP IDs)
SNPpos <- Tumor_LogR[,1:2]
Tumor_LogR <- Tumor_LogR[,3:114]
Blood_LogR <- Blood_LogR[,colnames(Tumor_LogR)]
Tumor_BAllele <- Tumor_BAllele[,colnames(Tumor_LogR)]
Blood_BAllele <- Blood_BAllele[,colnames(Tumor_LogR)]

# Divide the genome into 40 chromosome segments (chromosome arms, or when the centromere is small: whole chromosomes)
chr1p = 1:5213
chr1q = 5214:9819
chr2p = 9820:13130
chr2q = 13131:18521
chr3p = 18522:21782
chr3q = 21783:25728
chr4p = 25729:27475
chr4q = 27476:31728
chr5p = 31729:33220
chr5q = 33221:38057
chr6p = 38058:40940
chr6q = 40941:44636
chr7p = 44637:46798
chr7q = 46799:50217
chr8p = 50218:51941
chr8q = 51942:55108
chr9p = 55109:56619
chr9q = 56620:59588
chr10p = 59589:61027
chr10q = 61028:64828
chr11p = 64829:67122
chr11q = 67123:70756
chr12p = 70757:72305
chr12q = 72306:76221
chr13q = 76222:79314
chr14q = 79315:82734
chr15q = 82735:86041
chr16p = 86042:87543
chr16q = 87544:89429
chr17 = 89430:93508
chr18p = 93509:94073
chr18q = 94074:96078
chr19p = 96079:97630
chr19q = 97631:99598
chr20p = 99599:100822
chr20q = 100823:102605
chr21q = 102606:103986
chr22q = 103987:105872
chrXp = 105873:107264
chrXq = 107265:109302

chr = list(chr1p,chr1q,chr2p,chr2q,chr3p,chr3q,chr4p,chr4q,chr5p,chr5q,chr6p,chr6q,chr7p,chr7q,chr8p,chr8q,chr9p,chr9q,chr10p,chr10q,chr11p,chr11q,chr12p,chr12q,chr13q,chr14q,chr15q,chr16p,chr16q,chr17,chr18p,chr18q,chr19p,chr19q,chr20p,chr20q,chr21q,chr22q,chrXp,chrXq);


# make input files for the ASPCF algorithm
for (sample in 1:112) {
  tbsam = Tumor_BAllele[,sample]
  blsam = Blood_BAllele[,sample]
  Select_nonNAs <- blsam > 0.3 & blsam < 0.7 & !is.na(blsam) & !is.na(tbsam)
  tbsam[!Select_nonNAs] = -1
  for (chrke in 1:length(chr)) {
    lr = Tumor_LogR[chr[[chrke]],sample]
    baf = tbsam[chr[[chrke]]]
    lrNAs = is.na(lr)
    lr = lr[!lrNAs]
    baf = baf[!lrNAs]
    write.table(lr, paste(paste(paste("LRR",sample, sep=""),chrke,sep="_"),".txt", sep = ""), sep = "\t", row.names=F, col.names=F) 
    write.table(baf, paste(paste(paste("BAF",sample, sep=""),chrke,sep="_"),".txt", sep = ""), sep = "\t", row.names=F, col.names=F)
  }
}





# read output of ASPCF algorithm
Tumor_LogR_PCFed = matrix(nrow = dim(Tumor_LogR)[1], ncol = dim(Tumor_LogR)[2])
rownames(Tumor_LogR_PCFed) = rownames(Tumor_LogR)
Tumor_BAllele_PCFed = list();
for (sample in 1:112) {
  print(sample)
  lr = Tumor_LogR[,sample]
  lrNAs = is.na(lr)
  sampleData = matrix(nrow = length(na.omit(lr)), ncol = 2)
  rownames(sampleData) = rownames(Tumor_LogR)[!lrNAs]
  for (chrke in 1:length(chr)) {
    snps = rownames(Tumor_LogR)[chr[[chrke]]]
    lrb = Tumor_LogR[chr[[chrke]],sample]
    lrbNAs = is.na(lrb)
    snps = snps[!lrbNAs]
    lrchr = read.table(paste(paste(paste("LRR_PCF",sample, sep=""),chrke,sep="_"),".txt", sep = ""))[,1]
    bafchr = read.table(paste(paste(paste("BAF_PCF",sample, sep=""),chrke,sep="_"),".txt", sep = ""))[,1]
    sampleData[snps,1] = lrchr
    sampleData[snps,2] = bafchr
  }
  Tumor_LogR_PCFed[rownames(sampleData),sample] = sampleData[,1]
  tbsam = Tumor_BAllele[,sample]
  blsam = Blood_BAllele[,sample]
  Select_nonNAs <- blsam > 0.3 & blsam < 0.7 & !is.na(blsam) & !is.na(tbsam)
  rn = intersect(rownames(Tumor_BAllele)[Select_nonNAs],rownames(sampleData))
  Tumor_BAllele_PCFed[[sample]] = matrix(nrow = length(rn), ncol = 1)
  rownames(Tumor_BAllele_PCFed[[sample]]) = rn
  Tumor_BAllele_PCFed[[sample]][rn,1] = 1-sampleData[rn,2]
}



# Divide the genome into chromosomes (for visualization)
c1 = 1:9819
c2 = 9820:18521
c3 = 18522:25728
c4 = 25729:31728
c5 = 31729:38057
c6 = 38058:44636
c7 = 44637:50217
c8 = 50218:55108
c9 = 55109:59588
c10 = 59589:64828
c11 = 64829:70756
c12 = 70757:76221
c13 = 76222:79314
c14 = 79315:82734
c15 = 82735:86041
c16 = 86042:89429
c17 = 89430:93508
c18 = 93509:96078
c19 = 96079:99598
c20 = 99599:102605
c21 = 102606:103986
c22 = 103987:105872
cX = 105873:109302

ch = list(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,cX)


# plot the LRR and BAF (with probes homozygous in the germline removed) data (for tumor), before and after PCF
for (arraynr in 1:112) {
  png(filename = paste(paste("resPCF",arraynr,sep=""),".png",sep=""), width = 1500, height = 900)
  par(mfrow=c(2,1), mar = c(1,2,0.2,0.2), xaxt = 'n')
  r = Tumor_LogR_PCFed[rownames(Tumor_BAllele_PCFed[[arraynr]]),arraynr]
  beta = Tumor_BAllele_PCFed[[arraynr]][,]
  plot(c(1,length(r)), c(-1,1), type = "n", xlab = "", ylab = "")
  points(Tumor_LogR[rownames(Tumor_BAllele_PCFed[[arraynr]]),arraynr], col = "red")
  points(r,col="green")
  abline(v=0.5,lty=1,col="lightgrey")
  chrk_tot_len = 0
  for (i in 1:length(ch)) {
    chrk = ch[[i]];
    chrk_hetero = intersect(rownames(Tumor_LogR_PCFed)[chrk],rownames(Tumor_BAllele_PCFed[[arraynr]]))
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len+0.5;
    abline(v=vpos,lty=1,col="lightgrey")
  }
  plot(c(1,length(beta)), c(0,1), type = "n", xlab = "B Allele Frequency", ylab = "")
  points(Tumor_BAllele[rownames(Tumor_BAllele_PCFed[[arraynr]]),arraynr], col = "red")
  points(beta, col = "green")
  points(1-beta, col = "green")
  abline(v=0.5,lty=1,col="lightgrey")
  chrk_tot_len = 0
  for (i in 1:length(ch)) {
    chrk = ch[[i]];
    chrk_hetero = intersect(rownames(Tumor_LogR_PCFed)[chrk],rownames(Tumor_BAllele_PCFed[[arraynr]]))
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len+0.5;
    abline(v=vpos,lty=1,col="lightgrey")
  }
  dev.off()
}


# run the ASCAT algorithm for all arrays
source("ASCAT.R")
goodarrays=NULL
res = vector("list",112)
for (arraynr in 1:112) {
  print(arraynr)
  lrr=Tumor_LogR[,arraynr]
  names(lrr)=rownames(Tumor_LogR)
  baf=Tumor_BAllele[,arraynr]
  names(baf)=rownames(Tumor_BAllele)
  res[[arraynr]] = runASCAT(lrr,baf,Tumor_LogR_PCFed[,arraynr],Tumor_BAllele_PCFed[[arraynr]][,],ch,paste(paste("dist",arraynr,sep=""),".png",sep=""),paste(paste("profiles",arraynr,sep=""),".png",sep=""))
  if(!is.na(res[[arraynr]]$rho)) {
    goodarrays[length(goodarrays)+1] = arraynr
  }
}


# put all results into two matrices: n1 contains the copy number of the A allele, n2 contains the copy number of the B allele
n1 = matrix(nrow = dim(Tumor_LogR)[1], ncol = length(goodarrays))
n2 = matrix(nrow = dim(Tumor_LogR)[1], ncol = length(goodarrays))
rownames(n1) = rownames(Tumor_LogR)
rownames(n2) = rownames(Tumor_LogR)
colnames(n1) = colnames(Tumor_LogR)[goodarrays]
colnames(n2) = colnames(Tumor_LogR)[goodarrays]
for (i in 1:length(goodarrays)) {
  n1[,i] = res[[goodarrays[i]]]$nA
  n2[,i] = res[[goodarrays[i]]]$nB
}



