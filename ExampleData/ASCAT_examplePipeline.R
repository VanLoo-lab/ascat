#simplest form to run ASCAT, with matched normal data available, without GC wave correction and all samples female:

library(ASCAT)
ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt","Germline_LogR.txt","Germline_BAF.txt")
ascat.plotRawData(ascat.bc)
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)


#run ASCAT with additional GC correction (assuming AffySNP6 data, for other platforms please check our website)

library(ASCAT)
ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt","Germline_LogR.txt","Germline_BAF.txt")
ascat.bc <- ascat.GCcorrect(ascat.bc, "GC_AffySNP6_102015.txt")
ascat.plotRawData(ascat.bc)
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)

#run ASCAT without matched normals (assuming AffySNP6 data, for other platforms please check our website)

library(ASCAT)
ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt")
ascat.plotRawData(ascat.bc)
gg<-ascat.predictGermlineGenotypes(ascat.bc, platform = "AffySNP6")
ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=gg)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)

#run ASCAT with multi-sample segmentation (when shared breakpoints are expected)

library(ASCAT)
ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt","Germline_LogR.txt","Germline_BAF.txt")
ascat.plotRawData(ascat.bc)
ascat.bc = ascat.asmultipcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)
