# read raw methylation data

library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

load("data/target.RData")


#read Raw data
loadData <- function(target, out_dir){
  # read raw data
  rgSet <- read.metharray.exp(targets = target)
  
  # change sample name from sentrix information to real sample name (patient id - time)
  old_name <- sampleNames(rgSet)
  if(sum(old_name != paste(target$Sentrix_ID, target$Sentrix_Position, sep="_"))!=0){
    print("please check Sentrix id and Sentrix position in target file")
    return()
  }
  sampleNames(rgSet) <- target$Sample_Name
  
  detP <- detectionP(rgSet)
  pdf(file.path(out_dir, "detection.pdf"))
  barplot(colMeans(detP))
  dev.off()
  
  qcReport(rgSet, sampNames=target$Sample_Name, sampGroups=target$Sample_Group, 
           pdf = file.path(out_dir, "qcReport.pdf"))
  
  
  MSet <- preprocessRaw(rgSet)
  RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
  GRset <- mapToGenome(RSet)
  
  predictedSex <- getSex(GRset)
  
  rlt <- list(
    rgSet = rgSet,
    predictedSex = predictedSex,
    pvProb = detP
  )
  
  saveRDS(rlt, file.path(out_dir, "rawData.rds"))
}

loadData(target, "rlt/QC/")
