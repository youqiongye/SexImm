##R version 3.5, platform x86_64-redhat-linux-gnu
options(stringsAsFactors = F)
library(dummies)
source("/extraspace/yye1/analysis/2019SexImm/code/CodeSummary/cal.R")
##setup new folder
folder <- "3.4ImmunePopulations"
setwd("/extraspace/yye1/analysis/2019SexImm/3.PSM/")
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- paste("/extraspace/yye1/analysis/2019SexImm/3.PSM/",folder,sep="")
setwd(scripts.dir)
##The cancer types we investigated
cancerNames <- c("ACC", "BLCA", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG",
                 "LIHC", "LUAD", "LUSC", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD", "THCA", "THYM", "UVM")

###Loading the pan cancer immune cell populations abundance 
ImmuneSubsScore <- readr::read_rds(path="/extraspace/yye1/analysis/2019SexImm/code/CodeSummary/Data/PanCancerImmunePopulationAbundance.rds.gz")
sum.ImmuneAll <- data.frame()
for(cancer in cancerNames){
  data <- read.csv(paste("/extraspace/yye1/analysis/2019SexImm/ClinicalData/ProcessedClinicalData/",cancer,"_TumorPurity_ClinicalData.csv",sep = ""),header=T,stringsAsFactors=F)
  data$age_at_initial_pathologic_diagnosis <- as.numeric(data$age_at_initial_pathologic_diagnosis)
  data <- unique(data) 
  rownames(data) <- data[,1]
  data <- data[,-1]
  
  analysis <- "gender"
  if(analysis=="gender"){
    # convert female and male to numeric 1,0 to suppress the warning message in lm
    data$gender <- ifelse(data$gender=="FEMALE",1,0)
    colnames(data)[which(colnames(data)=="gender")] <- "Z"
  }
  
  # convert to dummy
  
  dummy.feature <- setdiff(colnames(data),c("Z","age_at_initial_pathologic_diagnosis","Purity"))#,"pathologic_stage"))
  if(length(dummy.feature)>0){
    data.dum <- dummy.data.frame(data, names=dummy.feature)
    dummy.list <- attr(data.dum,"dummies")
    rm.col <- c()
    for (i in 1:length(dummy.list))
    {
      rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
    }
    data.dum <- data.dum[,-rm.col]
    data.dum$X0 <- rep(1, nrow(data.dum))
    #form <- as.formula("Z~.") # should exclude X0
    exclude.col <- match(c("Z","X0"), colnames(data.dum))
    colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
    form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
  }else{
    data.dum <- data
    data.dum$X0 <- rep(1, nrow(data.dum))
    #form <- as.formula("Z~.") # should exclude X0
    exclude.col <- match(c("Z","X0"), colnames(data.dum))
    colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
    form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
  }
  
  
  # perform calculation
  
  
  # ImmuneCells
  EachImmune <- ImmuneSubsScore[ImmuneSubsScore$cancer_types==cancer,]$IS_GSVA[[1]]
  EachImmune <- EachImmune[substr(EachImmune$barcode,14,15) %in% c("01"),] %>% dplyr::mutate(barcode=substr(barcode,1,12))
  EachImmune <- EachImmune[!duplicated( EachImmune$barcode),]
  EachImmune.pri <- EachImmune[2:ncol(EachImmune)]
  colnames(EachImmune.pri) <- colnames(EachImmune)[2:ncol(EachImmune)]
  rownames(EachImmune.pri) <- EachImmune$barcode
  EachImmune.pri <- rm.zero.col(EachImmune.pri)
  
  EachImmune.result <- weight.test(data.dum, form, EachImmune.pri, is.continuous=TRUE,weight=ifelse(analysis=="gender","ATT","MW"),mirror.plot=FALSE, cancer, data.type= "EachImmune", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=FALSE)
  
  sum.Immune <- summarize.p(EachImmune.pri, EachImmune.result, print=TRUE)
  summarize.p(EachImmune.pri, EachImmune.result, print=TRUE, cutoff=0.05)
  write.summary(sum.Immune, cancer, analysis,"Immune")
  write.result(EachImmune.result, cancer, analysis,"Immune")
  save(EachImmune.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))
  if(length(which(EachImmune.result$pvalue < 0.05)) > 0){
    sum.Immune <- data.frame(sum.Immune)
    sum.Immune$class <- rep(cancer,times=nrow(sum.Immune))
    if(nrow(sum.ImmuneAll) == 0){
      sum.ImmuneAll <- sum.Immune
    }else{
      sum.ImmuneAll <- rbind(sum.ImmuneAll,sum.Immune)
    }
    
  }
}

write.table(sum.ImmuneAll,file="ImmuneCellsAbundance.across.cancer.typesAll.txt",quote = F,sep="\t",row.names = F)
