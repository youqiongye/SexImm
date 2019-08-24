###The association between the effect of immune checkpoint blockade (ICB) and gender within each trial
###There 27 ICB clinical trials obtained from Conforti et al (PMID: 29778737) and Wallis et al (PMID: 30605213)
options(stringsAsFactors = F)
library(magrittr)
library(metawho) ##Deft approached to assessing the effect of gender as measured within each relevant trial.

##Load clinical trial data, Data deposited in https://github.com/youqiongye/SexImm/tree/master/Data/
JLdata <- readxl::read_xlsx("Overall survival of ICB Clinic trials.xlsx")
##Cancer abbreviation organize
JLdata$Cancer <- ifelse(JLdata$Cancer %in% "gastroeso-oesophageal juncton cancer","GOJC",JLdata$Cancer)
JLdata$Cancer <- factor(JLdata$Cancer,levels=rev(c("Melanoma","NSCLC","SCLC","ccRCC","HNSC","Urothelial","GOJC","Mesothelioma")))

##get low border and upper border 95% CI
JLdata[,c("MEN_CIL","MEN_CIH")] <- apply(data.frame(do.call(rbind,strsplit(as.character(JLdata$MEN_CI),"-"))),2,as.numeric)
JLdata[,c("WOMEN_CIL","WOMEN_CIH")] <- apply(data.frame(do.call(rbind,strsplit(as.character(JLdata$WOMEN_CI),"-"))),2,as.numeric)

##Efficacy tendency for female and male patients with ICB treatment
JLdata$DiffClass <- ifelse(JLdata$HRDif > 0,"MenBetter","FemaleBetter")
JLdata <- JLdata[rev(order(JLdata$Cancer,JLdata$DiffClass,JLdata$WOMEN)),]
JLdata$Author <- factor(JLdata$Author)
JLdata_F <- JLdata[,c("Author","Cancer","WOMEN","WOMEN_CIL","WOMEN_CIH")] %>% set_colnames(c(c("Author","Cancer","HR","CIL","CIH"))) %>% 
  dplyr::mutate(gender = rep("Female",nrow(.)))
JLdata_M <- JLdata[,c("Author","Cancer","MEN","MEN_CIL","MEN_CIH")] %>% set_colnames(c(c("Author","Cancer","HR","CIL","CIH"))) %>% 
  dplyr::mutate(gender = rep("Male",nrow(.)))
JLdata.m <- rbind(JLdata_F,JLdata_M)
JLdata.m$Author <- factor(JLdata.m$Author,levels=JLdata$Author)
JLdata_PatientNum <- reshape2::melt(JLdata,id.vars="Author",measure.vars=c("WOMEN_Num","MEN_Num")) %>%
  set_colnames(c("Author","gender","value")) %>%
  dplyr::mutate(gender=ifelse(gender=="WOMEN_Num","Female","Male"))
###Forest plot of association between overall survival and immunotherapy (IO) and standard of care (SOC) stratified by patient sex
pdf("/extraspace/yye1/analysis/2019SexImm/ClinicalData/JAMA LancetOncologyData.pdf",width=8,height =4)
ggplot(JLdata.m,aes(y=HR,x= gender ,color = gender))+
  geom_segment(data= JLdata.m,aes(x=gender,y=CIL,xend=gender,yend=CIH),color="black",size=0.8,guide=F)+
  geom_point(size=3)+
  geom_text(data = JLdata_PatientNum[JLdata_PatientNum$gender=="Male",],aes(y=1.5,x= gender ,color = gender,label=paste("M (n = ",value,")",sep="")),angle = 90,vjust=0.5,hjust=0.5)+
  geom_text(data = JLdata_PatientNum[JLdata_PatientNum$gender=="Female",],aes(y=1.5,x= gender ,color = gender,label=paste("F (n = ",value,")",sep="")),angle = 90,vjust=0.5,hjust=0.5)+
  scale_color_manual(limits=c("Female","Male"),values=c("red","blue"),guide=F)+
  facet_grid(~Author)+
  theme(panel.background=element_rect(fill="gainsboro"),
        panel.grid=element_blank(),panel.spacing =  unit(0.1, "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),axis.line.y = element_line(color="black"),
        strip.text.x = element_text(angle=90,hjust = 0),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key = element_rect(fill="white",colour = "black"))
 dev.off()

##Heatmap to indicate the clinical trials invovled in which meta-analysis paper.
JLdata_paper <- reshape2::melt(JLdata,id.vars="Author",measure.vars=c("JAMA oncology","Lancet oncology"))
pdf("Clinical trials Resource.pdf",width=8,height =4)
ggplot(JLdata_paper,aes(y=variable,x= Author,fill = value)) +
  geom_tile(color="gray",size=1)+
  scale_x_discrete(limit=JLdata$Author)+
  coord_fixed()+
  scale_fill_manual(limit=c("Yes","No"),values=c("darkgreen","white"),guide=F)+
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.text.x=element_blank(),axis.title = element_blank(),axis.ticks = element_blank())
dev.off()
##

##Use of a deft approach to analyse and present interactions in meta-analysis, illustrating how the effect of ICB might vary by gender.
##### specify hazard ratios (hr)
hr <- paste(apply(JLdata[,c("MEN","WOMEN")],1,function(x){c(x[1],x[2])})) %>% as.numeric() 
### specify lower bound for hr confidence intervals
ci.lb <- paste(apply(JLdata[,c("MEN_CIL","WOMEN_CIL")],1,function(x){c(x[1],x[2])})) %>% as.numeric() 
### specify upper bound for hr confidence intervals
ci.ub <- paste(apply(JLdata[,c("MEN_CIH","WOMEN_CIH")],1,function(x){c(x[1],x[2])})) %>% as.numeric() 
### trials
trial <- rep(JLdata$Author,each=2)
### subgroups
subgroup <- rep(c("Male","Female"),nrow(JLdata))
entry <- paste(trial, subgroup, sep = "-")
MetaData =
  data.frame(
    entry = entry,
    trial = trial,
    subgroup = subgroup,
    hr = hr,
    ci.lb = ci.lb,
    ci.ub = ci.ub,
    stringsAsFactors = FALSE
  )

MetaData_pre <- deft_prepare(MetaData)
MetaData_pre$trial <- factor(MetaData_pre$trial,levels=JLdata$Author)
# The 'Male' is the reference
res = deft_do(MetaData_pre, group_level = c("Male", "Female"),method="REML")
# Extract the data after deft apporach, get hazard ratio (95% CI) between male and female patients comparison.
MetaData_Benefit <- res$subgroup$data
MetaData_Benefit$gender <- ifelse(MetaData_Benefit$hr < 1, "Female","Male")
MetaData_Benefit$trial <- factor(MetaData_Benefit$trial, levels=JLdata$Author)
MetaData_Benefit[,c("hr","ci.lb","ci.ub")] <- apply(MetaData_Benefit[,c("hr","ci.lb","ci.ub")],2,function(x)signif(x,digits = 2))

pdf("The efficacy of gender in ICB trials.pdf",width=8,height =3.5)
ggplot(MetaData_Benefit,aes(y=hr,x= 1 ,fill = gender))+
  geom_segment(data= MetaData_Benefit,aes(x=1,y=ci.lb,xend=1,yend=ci.ub),color="black",size=0.5,guide=F)+
  geom_point(size=3,shape=22)+
  scale_y_log10()+
  geom_text(data=MetaData_Benefit,aes(y=0.3,x=1,color=gender,label=paste(hr," (",ci.lb,"-",ci.ub,")",sep="")),angle = 90,vjust=0.5,hjust=0.5)+
  scale_fill_manual(limits=c("Female","Male"),values=c("red","blue"),guide=F)+
  scale_color_manual(limits=c("Female","Male"),values=c("red","blue"),guide=F)+
  
  facet_grid(~trial)+
  theme(panel.background=element_rect(fill="gainsboro"),
        panel.grid=element_blank(),panel.spacing =  unit(0.1, "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),axis.line.y = element_line(color="black"),
        strip.text.x = element_text(angle=90,hjust = 0),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()

