##The effect of confounders between female and male patients, 
##wilcox.test used for continous features, fihser.test used for discrete features.
options(stringsAsFactors = F)
library(magrittr)

setwd("./ProcessedClinicalData/")
###22 cancer types we investigated
cancerNames <- c("ACC", "BLCA", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG",
                 "LIHC", "LUAD", "LUSC", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD", "THCA", "THYM", "UVM")
###Two types of features
ContinousFeatures <- c("age_at_initial_pathologic_diagnosis","Purity")
DiscreteFeatures <- c("race","pathologic_stage","histological_type","tobacco_smoking_history")

##Statistic analysis
lapply(cancerNames,function(x){
  read.csv(paste(x,"_TumorPurity_ClinicalData.csv",sep="")) -> data
  Continous <- data[,c("gender",colnames(data)[colnames(data) %in% ContinousFeatures])]
  Continous[,2:ncol(Continous)] <- apply(Continous[,2:ncol(Continous)],2,function(x){ifelse(x=="#N/A",NA,x)})
  ContinousP <- data.frame()
  for(i in colnames(Continous)[-1]){
    ContinousP <- rbind(ContinousP,data.frame(Feature=i,pvalue= wilcox.test(as.numeric(Continous[Continous$gender=="FEMALE",i]),as.numeric(Continous[Continous$gender=="MALE",i]),na.action=na.omit)$p.val))
  }
  if(length(colnames(data)[colnames(data) %in% DiscreteFeatures])>0){
    Discrete <- data[,c("gender",colnames(data)[colnames(data) %in% DiscreteFeatures])]
    DiscreteP <- data.frame()
    for(i in colnames(Discrete)[-1]){
      DiscreteP <- rbind(DiscreteP,data.frame(Feature=i,pvalue= fisher.test(table(Discrete[,c("gender",i)]),simulate.p.value = TRUE)$p.value))
    }
    FeatureP <- rbind(ContinousP,DiscreteP)
    FeatureP$cancer_types <- rep(x,length(colnames(data)[colnames(data) %in% c(DiscreteFeatures,ContinousFeatures)]))
    return(FeatureP)
  }else{
    ContinousP$cancer_types <- rep(x,length(colnames(data)[colnames(data) %in% c(ContinousFeatures)]))
    return(ContinousP)
  }
}) %>% dplyr::bind_rows() -> FeaturesBias
#write.table(FeaturesBias,file="CounfounderAffectsBetweenGender.tab",quote = F,row.names = F,sep="\t")

FeaturesBias <- merge(FeaturesBias,data.frame(cancer_types=rep(cancerNames,times=6),Feature=rep(c(DiscreteFeatures,ContinousFeatures),each=22)),by=c("Feature","cancer_types"),all=T)
FeaturesBias$tmp <- FeaturesBias$pvalue
FeaturesBias[,"tmp"][is.na(FeaturesBias[,"tmp"])] <- 1
FeaturesBias[,"tmp"][FeaturesBias[,"tmp"] <= 0.05]  <- 0
FeaturesBias[,"tmp"][FeaturesBias[,"tmp"] > 0.05 & FeaturesBias[,"tmp"] <1]  <- 0.5
FeaturesBias$tmp <- factor(FeaturesBias$tmp)
pdf("./ProcessedClinicalData/Feature.bias.pdf",width=6,height = 2.5)#  order by immune cell abundance
ggplot(FeaturesBias[FeaturesBias$tmp %in% c(0,0.5),],aes(x=cancer_types,y=Feature))+
  geom_point(aes(color=tmp),size=4)+
  scale_color_manual(limit=c(0,0.5),values=c("#E41A1C","darkgray"),guide=F)+
  scale_y_discrete(limit=c(ContinousFeatures,DiscreteFeatures),labels=c("Age","Purity","Race","Stage","Histological type","Smoking"))+
  geom_point(data=FeaturesBias[FeaturesBias$tmp ==1,],aes(x=cancer_types,y=Feature),shape=4,size=4)+
  theme(panel.background=element_rect(colour="white",fill="white"),
        #panel.grid.major=element_line(color="black"),
        #panel.grid.minor=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(size=12,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.text.y=element_text(size=12,colour = "black"),
        axis.ticks=element_blank(),legend.text=element_text(size=14),
        axis.line=element_blank())+
  geom_tile(data=FeaturesBias,aes(x=cancer_types,y=Feature),fill=NA,color="black")
dev.off()
