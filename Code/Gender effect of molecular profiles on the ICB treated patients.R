##Gender effect of molecular biomarkers reported in immunotherapy datasets with molecular profiling for individual patients
##1. Data Set1 Miao et al.(PMID: 29301960). Molecular features: tumor mutation burden (TMB) and PBRM1 mutation. 
##2. Data Set2 Samstein et al. (PMID: 30643254). Molecular feature: TMB.
##3. Data Set3 Miao et al. (PMID: 30150660). Molecular features: TMB, neoantigen load, and APOBEC
##4. Data Set4 Cristescu et al. (PMID: 30309915). Molecular features: TMB and T cell–inflamed gene expression profile (GEP).
##5. Data Set5 Hugo et al. (PMID: 26997480). Molecular features: TMB and BRCA2 mutation
##6. Data Set6 Van Allen et al. (PMID: 26359337). Molecular features: TMB, neoantigen load, cytolytic activity (CYT), CTLA-4, and PD-L2
##7. Data Set7 Hellmann et al. (PMID: 29657128). Molecular features: TMB and PD-L1

############Package load, function define ###############
my.wilcox.test <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj[c("p.value")])
}
require("biomaRt")
library("magrittr")

############Process each dataset independently and merge result together to plot the figure.

##1. Data Set1 Miao et al.(PMID: 29301960). Molecular features: tumor mutation burden (TMB) and PBRM1 mutation. 
##Load TMB data
##
ccRCC_ClinicData <- read.csv("/extraspace/yye1/analysis/2019SexImm/5.ICBdata_FMdiff/5.1ccRCC_antiPD1/Raw_Data/ClinicalData_ccRCC_antiPD1.csv",stringsAsFactors = F)
DS1_ccRCC_TMB <- read.csv("DS1_PMID29301960_ccRCC_antiPD1.csv",stringsAsFactors = F)
DS1_ccRCC_TMB <- merge(DS1_ccRCC_TMB,ccRCC_ClinicData[,c("patient_id","sex")],by="patient_id")
DS1_ccRCC_TMBDiff <- data.frame(DataResource = "Science_ccRCC_antiPD1",Cancer.Type="ccRCC",
                            Female_TMB=median(log2(DS1_ccRCC_TMB[DS1_ccRCC_TMB$sex=="FEMALE",]$all_muts)),
                            Male_TMB= median(log2(DS1_ccRCC_TMB[DS1_ccRCC_TMB$sex=="MALE",]$all_muts)),
                            Diff = median(log2(DS1_ccRCC_TMB[DS1_ccRCC_TMB$sex=="FEMALE",]$all_muts))-median(log2(DS1_ccRCC_TMB[DS1_ccRCC_TMB$sex=="MALE",]$all_muts)),
                            Pval = unlist(my.wilcox.test(log2(DS1_ccRCC_TMB[DS1_ccRCC_TMB$sex=="FEMALE",]$all_muts),log2(DS1_ccRCC_TMB[DS1_ccRCC_TMB$sex=="MALE",]$all_muts))))

write.table(ccRCC_TMBDiff,file="/extraspace/yye1/analysis/2019SexImm/5.ICBdata_FMdiff/5.1ccRCC_antiPD1/ccRCC_TMBDiffSum.tab",quote = F,row.names=F,sep="\t")

pdf("/extraspace/yye1/analysis/2019SexImm/5.ICBdata_FMdiff/5.1ccRCC_antiPD1/ccRCC_TMBDiff between female and male.pdf",width = 3,height = 4)
ggplot(ccRCC_TMB,aes(x=sex,y=log2(all_muts+1),color=sex))+
  geom_boxplot(width=0.6,fill=NA,outlier.shape = NA)+
  ggbeeswarm::geom_quasirandom(alpha=0.5,width = 0.3)+
  scale_x_discrete(limit=c("MALE","FEMALE"))+ ylab("TMB") +
  scale_color_manual(limit= c("MALE","FEMALE"),values=c("blue","red"),name="")+
  theme(panel.background=element_rect(colour=NA,fill="lightgray",size=0.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),panel.spacing = unit(0.1,"line"),
        axis.text.x=element_text(color="black"),axis.title.y=element_text(size=12,color="black"), axis.title.x=element_blank(),
        axis.text.y=element_text(size=10,color="black"), axis.line=element_line(),axis.ticks.x = element_blank(),
        axis.ticks.length = unit(.15, "cm"),strip.text=element_text(size=10),
        strip.background = element_rect(fill=NA),legend.position = "bottom")+
  geom_text(data = ccRCC_TMBDiff,aes(x=.5,y=8,label=paste("Diff = ",signif(Diff,digits = 2),"\np = ",signif(Pval,digits = 2),sep="")),color="black",size=4,hjust=0)
dev.off()
##2. Data Set2 Samstein et al. (PMID: 30643254). Molecular feature: TMB.
##3. Data Set3 Miao et al. (PMID: 30150660). Molecular features: TMB, neoantigen load, and APOBEC
##4. Data Set4 Cristescu et al. (PMID: 30309915). Molecular features: TMB and T cell–inflamed gene expression profile (GEP).
##5. Data Set5 Hugo et al. (PMID: 26997480). Molecular features: TMB and BRCA2 mutation
##6. Data Set6 Van Allen et al. (PMID: 26359337). Molecular features: TMB, neoantigen load, cytolytic activity (CYT), CTLA-4, and PD-L2
##7. Data Set7 Hellmann et al. (PMID: 29657128). Molecular features: TMB and PD-L1