#############################################
# First, CEL files are preprocessed using PennCNV-Affy (giving LogR and BAF), see http://penncnv.openbioinformatics.org/en/latest/user-guide/affy/ only using steps 1.1, 1.2 and 1.4
# step 1.1:
# apt-probeset-genotype -c GenomeWideSNP_6.cdf -a birdseed --read-models-birdseed GenomeWideSNP_6.birdseed.models --special-snps GenomeWideSNP_6.specialSNPs --out-dir apt --cel-files CELfiles.txt
# step 1.2:
# apt-probeset-summarize --cdf-file GenomeWideSNP_6.cdf --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch hapmap.quant-norm.normalization-target.txt --out-dir apt --cel-files CELfiles.txt
# step 1.4:
# normalize_affy_geno_cluster.pl gw6.genocluster quant-norm.pm-only.med-polish.expr.summary.txt -locfile affygw6.hg19.pfb -out lrr_baf.txt
#############################################

lrrbaf = read.table("lrr_baf.txt", header = T, sep = "\t", row.names=1)

SNPpos = read.table("SNPpos.txt",header=T,sep="\t",row.names=1)

sample = sub(".normal.CEL.Log.R.Ratio","",colnames(lrrbaf)[3])

Tumor_LogR = lrrbaf[rownames(SNPpos),5,drop=F]
colnames(Tumor_LogR) = sample

Tumor_BAF = lrrbaf[rownames(SNPpos),6,drop=F]
colnames(Tumor_BAF) = sample

Normal_LogR = lrrbaf[rownames(SNPpos),3,drop=F]
colnames(Normal_LogR) = sample

Normal_BAF = lrrbaf[rownames(SNPpos),4,drop=F]
colnames(Normal_BAF) = sample

#replace 2's by NA
Tumor_BAF[Tumor_BAF==2]=NA
Normal_BAF[Normal_BAF==2]=NA

# Tumor_LogR: correct difference between copy number only probes and other probes
CNprobes = substring(rownames(SNPpos),1,2)=="CN"

Tumor_LogR[CNprobes,1] = Tumor_LogR[CNprobes,1]-mean(Tumor_LogR[CNprobes,1],na.rm=T)
Tumor_LogR[!CNprobes,1] = Tumor_LogR[!CNprobes,1]-mean(Tumor_LogR[!CNprobes,1],na.rm=T)

Normal_LogR[CNprobes,1] = Normal_LogR[CNprobes,1]-mean(Normal_LogR[CNprobes,1],na.rm=T)
Normal_LogR[!CNprobes,1] = Normal_LogR[!CNprobes,1]-mean(Normal_LogR[!CNprobes,1],na.rm=T)

# limit the number of digits:
Tumor_LogR = round(Tumor_LogR,4)
Normal_LogR = round(Normal_LogR,4)

write.table(cbind(SNPpos,Tumor_BAF),paste(sample, ".tumor.BAF.txt", sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
write.table(cbind(SNPpos,Normal_BAF),paste(sample, ".normal.BAF.txt", sep=""),sep="\t",row.names=T,col.names=NA,quote=F)

write.table(cbind(SNPpos,Tumor_LogR),paste(sample, ".tumor.LogR.txt", sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
write.table(cbind(SNPpos,Normal_LogR),paste(sample, ".normal.LogR.txt", sep=""),sep="\t",row.names=T,col.names=NA,quote=F)

#run ASCAT functions

library(ASCAT)
file.tumor.LogR <- dir(pattern="tumor.LogR")
file.tumor.BAF <- dir(pattern="tumor.BAF")
file.normal.LogR <- dir(pattern="normal.LogR")
file.normal.BAF <- dir(pattern="normal.BAF")

gender <- read.table("birdseed.report.txt", sep="\t", skip=66, header=T)
sex <- as.vector(gender[,"computed_gender"])
sex[sex == "female"] <- "XX"
sex[sex == "male"] <- "XY"
sex[sex == "unknown"] <- "XX"

samplename <- sub(".tumor.LogR.txt", "", file.tumor.LogR)

ascat.bc <- ascat.loadData(file.tumor.LogR, file.tumor.BAF, file.normal.LogR, file.normal.BAF, chrs=c(1:22, "X"), gender=sex)

#GC correction for SNP6 data
ascat.bc <- ascat.GCcorrect(ascat.bc, "GC_AffySNP6_102015.txt")

ascat.plotRawData(ascat.bc)

ascat.bc <- ascat.aspcf(ascat.bc)

ascat.plotSegmentedData(ascat.bc)

ascat.output <- ascat.runAscat(ascat.bc)

#save ASCAT results

write.table(ascat.output$segments, file=paste(samplename,".segments.txt",sep=""), sep="\t", quote=F, row.names=F)
write.table(ascat.output$aberrantcellfraction, file=paste(samplename,".acf.txt",sep=""), sep="\t", quote=F, row.names=F)
write.table(ascat.output$ploidy, file=paste(samplename,".ploidy.txt",sep=""), sep="\t", quote=F, row.names=F)

save.image(paste(samplename,".RData",sep=""))