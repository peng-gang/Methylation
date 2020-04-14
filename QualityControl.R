# Quality Control

library(minfi)
library(ChAMP)
library(ggplot2)
library(flashpcaR)
library(ggsci)
library(wateRmelon)
library(FlowSorted.Blood.EPIC)

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# load data
load("data/target.RData")

raw <- readRDS("rlt/QC/rawData.rds")

# number of low quality probes for each sample
# max is 809 out of 866K
sort(colSums(raw$pvProb>0.01))

probFaile <- rowSums(raw$pvProb>0.01)

idxLowQualityLoci <- probFaile > nrow(target) * 0.05
# 1737 probes with low quality
sum(idxLowQualityLoci)

lowQualityLoci <- rownames(raw$pvProb)[idxLowQualityLoci]

# there might be warning from minfi package if we remove low quality probe here.
# we remove low quality probe after preprocessQuantile


# sex match
sum(raw$predictedSex[,3] == ifelse(target$Sex=="Female", "F", "M"))


# quantile normalization
GR <- preprocessQuantile(
  raw$rgSet, fixOutliers = TRUE, 
  removeBadSamples = TRUE,
  badSampleCutoff = 10.5,
  quantileNormalize = TRUE,
  stratified = TRUE,
  mergeManifest = TRUE
)

# remove SNP related CpG sites
GR <- dropLociWithSnps(GR, snps=c("SBE","CpG"), maf=0)

# remove low quality, multiple hit, chrY probes (all female)

## load multi.hit data from package "ChAMP"
data(multi.hit)

probeAll <- featureNames(GR)

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

idxLowQualityProbe <- probeAll %in% lowQualityLoci
idxMultiProbe <- probeAll %in% as.character(multi.hit$TargetID)
idxChrYProbe <- probeAll %in% annEPIC$Name[annEPIC$chr %in% "chrY"]

idxRmAll <- idxLowQualityProbe | idxMultiProbe | idxChrYProbe

sum(idxRmAll)

GR <- GR[!idxRmAll,]

save(GR, file = "data/GR.RData")


# quality check
M <- getM(GR)
Beta <- getBeta(GR)

varBeta <- apply(Beta, 1, var)
idxVarBeta <- order(varBeta, decreasing = TRUE)

hist(varBeta)

#MDS
#BF10073-3, BF10049-3, BF10039-2 in race 1, others in race 2
limma::plotMDS(M, top = 1000, gene.selection = "common",
               col = ifelse(target$Sample_Group == "affected", "red", "blue"))

limma::plotMDS(Beta, top = 1000, gene.selection = "common",
               col = ifelse(target$Sample_Group == "affected", "red", "blue"))


limma::plotMDS(M, top = 10000, gene.selection = "common",
               col = ifelse(target$Sample_Group == "affected", "red", "blue"))

limma::plotMDS(Beta, top = 10000, gene.selection = "common",
               col = ifelse(target$Sample_Group == "affected", "red", "blue"))



#PCA
plotPCA <- function(x, group, scale = FALSE, top=0, varAll = NULL){
  if(top <=0 | top >= nrow(x)){
    xScale <- scale(t(x), center = TRUE, scale = scale)
    pcaRaw <- flashpca(xScale, ndim = 5, stand = "none", divisor = "none")
    
    dplot <- data.frame(
      x = pcaRaw$vectors[,1], y = pcaRaw$vectors[,2],
      group = factor(group)
    )
    
    gp <- ggplot(dplot) + 
      geom_point(aes(x=x, y=y, color=group)) + 
      labs(x="PC1", y = "PC2") + 
      scale_color_npg() + 
      theme_light() + 
      theme(text = element_text(size=20))
    return(gp)
    
  } else {
    if(is.null(varAll)){
      varAll <- apply(x, 1, var)
    }
    
    varOrder <- order(varAll, decreasing = TRUE)
    xScale <- scale(t(x[varOrder[1:top],]), center = TRUE, scale = scale)
    
    pcaRaw <- flashpca(xScale, ndim = 5, stand = "none", divisor = "none")
    
    dplot <- data.frame(
      x = pcaRaw$vectors[,1], y = pcaRaw$vectors[,2],
      group = factor(group)
    )
    
    gp <- ggplot(dplot) + 
      geom_point(aes(x=x, y=y, color=group)) + 
      labs(x="PC1", y = "PC2") + 
      scale_color_npg() + 
      theme_light() + 
      theme(text = element_text(size=20))
    return(gp)
  }
}



gp <- plotPCA(Beta, target$Sample_Group, scale = FALSE, top=1000, varBeta)
gp <- gp + theme(legend.position = c(0.8, 0.8), legend.title = element_blank())

pdf("rlt/PCATop500.pdf", width = 5, height = 5)
gp
dev.off()
#plotPCA(Beta, target$Sample_Group, scale = TRUE, top=1000, varBeta)

gp <- plotPCA(Beta, target$Sample_Group, scale = FALSE, top=0, varBeta)
gp <- gp + theme(legend.position = c(0.5, 0.85), legend.title = element_blank())
pdf("rlt/PCAAll.pdf", width = 5, height = 5)
gp
dev.off()

#plotPCA(Beta, target$Sample_Group, scale = TRUE, top=0, varBeta)



# cell type composition cannot be applied to this dataset 
# it is from saliva
# cellCounts <- estimateCellCounts2(
#   raw$rgSet, 
#   compositeCellType = "Blood",
#   processMethod = "preprocessNoob",
#   probeSelect = "IDOL",
#   cellTypes = c("CD8T", "CD4T", "NK", "Bcell",
#                 "Mono", "Neu"),
#   referencePlatform =
#     "IlluminaHumanMethylationEPIC",
#   referenceset = NULL,
#   IDOLOptimizedCpGs =IDOLOptimizedCpGs,
#   returnAll = TRUE
# )

predictedAge <- agep(Beta)

dplot <- data.frame(
  age = target$Age,
  predictAge = predictedAge[,1],
  group = target$Sample_Group,
  stringsAsFactors = FALSE
)

gp <- ggplot(dplot) + 
  geom_point(aes(x=age, y = predictAge, color = group)) + 
  geom_abline(slope = 1, intercept = 0, color = "grey") + 
  scale_x_continuous(limits = c(10, 70)) +
  scale_y_continuous(limits = c(10, 70)) +
  labs(x="Real Age", y = "Predicted Age") + 
  theme_light() + scale_color_npg() + 
  theme(legend.position = c(0.8, 0.2), legend.title = element_blank()) + 
  theme(text = element_text(size=16))

pdf("rlt/predictedAge.pdf", width = 5, height = 5)
gp
dev.off()


dplot <- data.frame(
  age = target$Age,
  group = target$Sample_Group,
  stringsAsFactors = FALSE
)

gp <- ggplot(dplot) +
  geom_boxplot(aes(x=group, y = age)) + 
  theme_classic() + 
  labs(x="", y = "Age") +
  theme(text = element_text(size=16))

pdf("rlt/ageDiff.pdf", width = 5, height = 5)
gp
dev.off()

