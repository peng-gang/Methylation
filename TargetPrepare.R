# prepare target data
library(readxl)

sampleInfo <- read_xlsx("data/methylation - cases - controls.xlsx")
experimentDesign <- read_xlsx("data/Raw/WenzhongLiu_mendelian_EPIC_030920.xlsx")

sum(experimentDesign$`Sample ID` == experimentDesign$ID)
sum(sampleInfo$ID ==  experimentDesign$ID)

target <- data.frame(
  Sample_Name = experimentDesign$ID,
  Sample_Well = experimentDesign$well,
  Sample_Plate = "Plate", # 24 samples are in one plate
  Sample_Group = sampleInfo$`Affected/unaffected`,
  Pool_ID = NA,
  Sentrix_ID = experimentDesign$`Chip Barcode`,
  Sentrix_Position = experimentDesign$`Sentrix Position`,
  Basename = file.path("data/Raw", 
                       paste(experimentDesign$`Chip Barcode`, experimentDesign$`Sentrix Position`, sep="_")),
  Sex = sampleInfo$Sex,
  Age = sampleInfo$Age,
  stringsAsFactors = FALSE
)


save(target, file = "data/target.RData")
