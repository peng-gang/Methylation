# compare between affected and unaffected

library(ggplot2)
library(reshape2)
library(ggsci)
library(minfi)
library(readxl)
library(stringr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(DMRcate)


load("data/GR.RData")
load("data/target.RData")

# methylation
M <- getM(GR)
Beta <- getBeta(GR)

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


# probe level
plotMethylation <- function(probe, group, file = NULL){
  idx <- which(rownames(Beta) == probe)
  
  dplot <- data.frame(
    x = group,
    y = Beta[idx,],
    stringsAsFactors = FALSE
  )
  
  gp <- ggplot(dplot) + 
    geom_boxplot(aes(x=x, y=y)) + 
    labs(x = "", y = "Beta", title = probe) + 
    theme_light() + theme(plot.title = element_text(hjust = 0.5))
  
  if(is.null(file)){
    gp
  } else {
    pdf(file, width = 5, height = 5)
    print(gp)
    dev.off()
  }
}

group <- target$Sample_Group

# unadjusted model
design <- model.matrix(~group)
fit <- limma::lmFit(M, design)
fit2 <- limma::eBayes(fit)

probes <- limma::topTable(fit2, coef = "groupunaffected", p.value = 0.05, number = 1000000)
probesAll <- limma::topTable(fit2, coef = "groupunaffected", number = 1000000)

sigCpGs <- rownames(probes)

idx <- match(sigCpGs, annEPIC$Name)
rlt <- data.frame(annEPIC[idx, c(1:4, 18:19, 22:24)], 
                  pvalue = probes$P.Value, fdr = probes$adj.P.Val,
                  stringsAsFactors = FALSE)

write.csv(rlt[1:1000,], file = "rlt/unadjust/probe/PVTop1000.csv",
          row.names = FALSE, quote = FALSE)

for(i in 1:10){
  plotMethylation(
    rownames(probes)[i], group, 
    file.path("rlt/unadjust/probe",
              paste0(rownames(probes)[i], ".pdf"))
  )
}


# DMR
DMRAnnotation <- cpg.annotate(
  object = M,
  what = "M",
  datatype = "array",
  arraytype = "EPIC",
  analysis.type = "differential",
  design = design, coef = "groupunaffected"
)

DMRs <- dmrcate(DMRAnnotation)

dmrProbes <- NULL
for(i in 1:nrow(DMRs$results)){
  coord <- DMRs$results[i,1]
  chr_choose <- strsplit(coord, ":")[[1]][1]
  st <- as.numeric(strsplit(strsplit(coord, ":")[[1]][2], "-")[[1]][1])
  ed <- as.numeric(strsplit(strsplit(coord, ":")[[1]][2], "-")[[1]][2])
  idx <- annEPIC$chr==chr_choose & annEPIC$pos >= st & annEPIC$pos <= ed
  ann_sel <-  annEPIC[idx,]
  idx_b <- rownames(Beta) %in% ann_sel$Name
  b_val_sel <- Beta[idx_b,]
  probe_sel <- rownames(b_val_sel)
  dmrProbes <- c(dmrProbes, paste(probe_sel, collapse =";"))
}

DMR_rlt <- data.frame(DMRs$results[,c(1, 2, 4)], dmrProbes)

write.csv(DMR_rlt, file = "rlt/unadjust/DMR/DMRs.csv",
          row.names = FALSE, quote = FALSE)


#GSEA

allCpGs <- rownames(Beta)

gseaGO <- missMethyl::gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs, 
  collection = "GO",
  array.type = "EPIC",
  plot.bias = TRUE
)

gseaRlt <- missMethyl::topGSA(gseaGO, number = 200)
write.csv(gseaRlt, file = "rlt/unadjust/GSEA/GOTop200.csv",
          row.names = FALSE, quote = TRUE)


gseaKEGG <- missMethyl::gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs, 
  collection = "KEGG",
  array.type = "EPIC",
  plot.bias = TRUE
)

gseaRlt <- missMethyl::topGSA(gseaKEGG, number = 200)
write.csv(gseaRlt, file = "rlt/unadjust/GSEA/KeggTop200.csv",
          row.names = FALSE, quote = TRUE)











# adjust with age
age <- target$Age

design <- model.matrix(~group + age)
fit <- limma::lmFit(M, design)
fit2 <- limma::eBayes(fit)

# no significant probes after multiple test correction
# take all probes with unadjust p-values < 0.05
probes <- limma::topTable(fit2, coef = "groupunaffected", number = 159348)
probesAll <- limma::topTable(fit2, coef = "groupunaffected", number = 1000000)

sigCpGs <- rownames(probes)

idx <- match(sigCpGs, annEPIC$Name)
rlt <- data.frame(annEPIC[idx, c(1:4, 18:19, 22:24)], 
                  pvalue = probes$P.Value, fdr = probes$adj.P.Val,
                  stringsAsFactors = FALSE)

write.csv(rlt[1:1000,], file = "rlt/adjust/probe/PVTop1000.csv",
          row.names = FALSE, quote = FALSE)

for(i in 1:10){
  plotMethylation(
    rownames(probes)[i], group, 
    file.path("rlt/adjust/probe",
              paste0(rownames(probes)[i], ".pdf"))
  )
}


# DMR
DMRAnnotation <- cpg.annotate(
  object = M,
  what = "M",
  datatype = "array",
  arraytype = "EPIC",
  analysis.type = "differential",
  design = design, coef = "groupunaffected"
)

DMRs <- dmrcate(DMRAnnotation, pcutoff = 0.05)

dmrProbes <- NULL
for(i in 1:nrow(DMRs$results)){
  coord <- DMRs$results[i,1]
  chr_choose <- strsplit(coord, ":")[[1]][1]
  st <- as.numeric(strsplit(strsplit(coord, ":")[[1]][2], "-")[[1]][1])
  ed <- as.numeric(strsplit(strsplit(coord, ":")[[1]][2], "-")[[1]][2])
  idx <- annEPIC$chr==chr_choose & annEPIC$pos >= st & annEPIC$pos <= ed
  ann_sel <-  annEPIC[idx,]
  idx_b <- rownames(Beta) %in% ann_sel$Name
  b_val_sel <- Beta[idx_b,]
  probe_sel <- rownames(b_val_sel)
  dmrProbes <- c(dmrProbes, paste(probe_sel, collapse =";"))
}

DMR_rlt <- data.frame(DMRs$results[,c(1, 2, 4)], dmrProbes)

write.csv(DMR_rlt, file = "rlt/adjust/DMR/DMRs.csv",
          row.names = FALSE, quote = FALSE)


#GSEA

allCpGs <- rownames(Beta)

gseaGO <- missMethyl::gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs, 
  collection = "GO",
  array.type = "EPIC",
  plot.bias = TRUE
)

gseaRlt <- missMethyl::topGSA(gseaGO, number = 200)
write.csv(gseaRlt, file = "rlt/adjust/GSEA/GOTop200.csv",
          row.names = FALSE, quote = TRUE)


gseaKEGG <- missMethyl::gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs, 
  collection = "KEGG",
  array.type = "EPIC",
  plot.bias = TRUE
)

gseaRlt <- missMethyl::topGSA(gseaKEGG, number = 200)
write.csv(gseaRlt, file = "rlt/adjust/GSEA/KeggTop200.csv",
          row.names = FALSE, quote = TRUE)




## compare adjusted and unadjusted results
rltUnadj <- read.csv("rlt/unadjust/probe/PVTop1000.csv")
rltAdj <- read.csv("rlt/adjust/probe/PVTop1000.csv")

sum(rltUnadj$Name %in% rltAdj$Name)
sum(rltAdj$Name %in% rltUnadj$Name)

