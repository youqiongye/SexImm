##R version 3.5, platform x86_64-redhat-linux-gnu
options(stringsAsFactors = F)
library(dummies)
#load propensity score algorithm function.
source("cal.R") ##Code can downloaded from https://github.com/youqiongye/SexImm/edit/master/Data/.
##setup new folder, deposit result in this folder
folder <- "PSM_Output"
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- folder
setwd(scripts.dir)
##22 cancer types we investigated
cancerNames <- c("ACC", "BLCA", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG",
                 "LIHC", "LUAD", "LUSC", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD", "THCA", "THYM", "UVM")

###Loading demo data: the pan cancer mRNA expression of immune checkpoint
##Data can be downloaded from https://github.com/youqiongye/SexImm/edit/master/Data/
ImmFeature <- readr::read_rds("InputDataForPSM.rds") 

analysis="gender" 
sum.ImmFeatureAll <- data.frame()

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
  Feature <- ImmFeature[ImmFeature$cancer_types==cancer,]$InAcMarkerExp[[1]]
  Feature <- Feature[substr(Feature$barcode,14,15) %in% c("01"),] %>% dplyr::mutate(barcode=substr(barcode,1,12))
  Feature <- Feature[!duplicated( Feature$barcode),]
  Feature.pri <- Feature[2:ncol(Feature)]
  colnames(Feature.pri) <- colnames(Feature)[2:ncol(Feature)]
  rownames(Feature.pri) <- Feature$barcode
  Feature.pri <- rm.zero.col(Feature.pri)
  Feature.pri <- apply(Feature.pri,2,function(x){log2(x+1)})
  folder <- paste(cancer,"_",analysis,sep="")
  if (!file.exists(folder)) { dir.create(folder) }
  
  Feature.result <- weight.test(data.dum, form, Feature.pri, is.continuous=TRUE,weight=ifelse(analysis=="gender","ATT","MW"),mirror.plot=FALSE, cancer, data.type= "Feature", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=FALSE)
  
  sum.Immune <- summarize.p(Feature.pri, Feature.result, print=TRUE)
  summarize.p(Feature.pri, Feature.result, print=TRUE, cutoff=0.05)
  write.summary(sum.Immune, cancer, analysis,"Immune")
  write.result(Feature.result, cancer, analysis,"Immune")
  save(Feature.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))
  if(length(which(Feature.result$pvalue < 0.05)) > 0){
    sum.Immune <- data.frame(sum.Immune)
    sum.Immune$class <- rep(cancer,times=nrow(sum.Immune))
    if(nrow(sum.ImmFeatureAll) == 0){
      sum.ImmFeatureAll <- sum.Immune
    }else{
      sum.ImmFeatureAll <- rbind(sum.ImmFeatureAll,sum.Immune)
    }
    
  }
}

write.table(sum.ImmFeatureAll,file="feature difference.across.cancer.typesAll.txt",quote = F,sep="\t",row.names = F)
