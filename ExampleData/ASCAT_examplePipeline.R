#simplest form to run ASCAT, with matched normal data available, without GC wave correction and all samples female:

library(ASCAT)
ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt","Germline_LogR.txt","Germline_BAF.txt")
ascat.plotRawData(ascat.bc)
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)
