library(dplyr)

#######################################################################
##############################Influential Genes Identification#########
#######################################################################

setwd("F:\\Green University\\Personal\\MSc in RUET\\MSc Thesis\\Thesis Work\\Journal\\Influential Genes and Survival")
clinical_data <- read.csv("coadread_tcga_pan_can_atlas_2018_clinical_data.csv",header = TRUE,stringsAsFactors = FALSE)

clinical_data <- dplyr::rename(clinical_data,Stage=Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)
early_stage <- subset(clinical_data,Stage=="STAGE I")
late_stage <- subset(clinical_data,Stage!="STAGE I" & Stage!="NA")

zscore_data <- read.csv("data_mrna_seq_v2_rsem_zscores_ref_all_samples.csv",header = FALSE,stringsAsFactors = FALSE)

zscore_data <- zscore_data[rowSums(is.na(zscore_data))<ncol(zscore_data)/5,]


group1 <- matrix(nrow = nrow(zscore_data), ncol =nrow(early_stage))
group2 <- matrix(nrow = nrow(zscore_data), ncol =nrow(late_stage))
group1[,1] <- zscore_data[,1]
group2[,1] <- zscore_data[,1]

k <- 2
l <- 2

for(i in 2:ncol(zscore_data))
{
 
  for(j in 1:nrow(early_stage))
    if(zscore_data[1,i]==early_stage[j,2])
    {
      group1 [,k] <- zscore_data[,i]
      k <- k+1
     break
    }
}

for(i in 2:ncol(zscore_data))
{
  
  for(j in 1:nrow(late_stage))
    if(zscore_data[1,i]==late_stage[j,2])
    {
      group2 [,l] <- zscore_data[,i]
      l <- l+1
      break
    }
}

group1 <- group1[,colSums(is.na(group1))<nrow(group1)]
group2 <- group2[,colSums(is.na(group2))<nrow(group2)]

#write.table(group1, file = "group1.csv",sep=",",  row.names=FALSE, col.names=FALSE)
#write.table(group2, file = "group2.csv",sep=",",  row.names=FALSE, col.names=FALSE)



alpha <- 0.05

cnt <- 0

names <- vector('character')
pvals <- vector('numeric')
LogFC <- vector('numeric')

for (i in 2:nrow(zscore_data)){
  
  vec1 = as.numeric(group1[i, 2:ncol(group1)])
  vec2 = as.numeric(group2[i, 2:ncol(group2)])
  
  FC <- mean(vec1)/mean(vec2)
  LFC <- log2(abs(FC))
  
  if(FC<0)
    LFC <- -LFC
  
  tt = kruskal.test(list(vec1, vec2))
  names <- c(names, zscore_data[i,1])
  pvals <- c(pvals, tt$p.value)
  LogFC <- c(LogFC,LFC)
}

padjusted = p.adjust(pvals, method = "bonferroni")

res = data.frame("Gene"=names, "p-value"=pvals, "Adjusted p-value"=padjusted, "LogFC"=LogFC)

kept = padjusted < 0.05


Key_genes = data.frame("Gene"=names[kept], "p-value"=pvals[kept], "Adjusted p-value"=padjusted[kept],"LogFC"=LogFC[kept])

cat(nrow(Key_genes))

#write.csv(res, file="kruskal_test_result_all.csv",row.names=FALSE)
write.csv(Key_genes, file="Key_genes.csv",row.names=FALSE)


#######################################################################
##############################Survival Analysis of the Key Genes#######
#######################################################################

library(stringr) 
library(plyr)


##############################Preprocessing Z-scores data#########
Key_genes<-read.csv("Key_genes.csv",header = TRUE,stringsAsFactors = FALSE)

Zscores<-read.csv("data_mrna_seq_v2_rsem_zscores_ref_all_samples.csv",header = TRUE,stringsAsFactors = FALSE)
Zscores2 <- Zscores[rowSums(is.na(Zscores))<ncol(Zscores)/5,]

Zscores2<-Zscores[(Zscores$Hugo_Symbol %in% Key_genes$Gene),]

Zscores2<-dplyr::rename(Zscores2,Patient_ID=Hugo_Symbol)
Zscores2<-t(Zscores2)


write.table(Zscores2, file = "data_mrna_seq_v2_rsem_zscores_ref_all_samples_2.csv",sep=",",  col.names=FALSE)



##############################working with clinical data#########

dato<-read.csv("coadread_tcga_pan_can_atlas_2018_clinical_data.csv",header = TRUE,stringsAsFactors = FALSE)


datos<-select(dato,Patient.ID,Race.Category,Diagnosis.Age,Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code,Censor.Status,Overall.Survival.in.Days)


#ranme the column
datos<-dplyr::rename(datos,Patient_ID=Patient.ID,Rectime=Overall.Survival.in.Days,Age=Diagnosis.Age,Cancer_stage=Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code,Race=Race.Category,Censor_status=Censor.Status)


datos$Race=as.character(lapply(datos$Race,function(x){gsub("^$","Others",x)}))
#replace "" in cancer_stage with NA
datos$Cancer_stage=as.character(lapply(datos$Cancer_stage,function(x){gsub("^$",NA,x)}))


for(i in 1:nrow(datos))
{
  datos[i,1] <- str_replace_all(datos[i,1],"-",".")## Replacing charecter - by .
}

count(datos,'Cancer_stage')
count(datos,'Race')



#################################################################
##############################working with mRNA data#########
#################################################################



#read rna expression data
Zscores2<-read.csv("data_mrna_seq_v2_rsem_zscores_ref_all_samples_2.csv",header=TRUE,stringsAsFactors = FALSE)


#function for labelling each expression value
altered_test<-function(x){
  
  if(typeof(x)=="character"){
    d=x
  }
  else{
    
    if (abs(x)>=2){
      d="Altered"
      
    }
    else{
      d="Normal"
      
    }
  }
  d
}

#apply the function over all the colulmn to convert altered unaltered

applyfunc<-function(df,f){
  ds<-matrix(0,nrow = nrow(df),ncol=ncol(df))
  colnames(ds)<-colnames(df)
  for (i in seq(1:ncol(df))){
    ds[,i]<-(sapply(Zscores2[,i],f))
  }
  ds<-as.data.frame(ds)
}
gene_status<-applyfunc(Zscores2,altered_test)

#remove the 01 from patient iD

remove_01<-function(x){
  x<-unlist(strsplit(x,split=""))
  x<-paste(x[0:(length(x)-3)],collapse = "")
  x
}

#gene_status$Patient_ID<-as.character(gene_status$Patient_ID)
gene_status$Patient_ID=unlist(lapply(gene_status$Patient_ID,remove_01))
#View(gene_status)

#########################Merge the tables###########################

#gene_status$Patient_ID=as.character(gene_status$Patient_ID)
combined<-datos%>%inner_join(gene_status)
#View(combined)  
#View(datos)
#View(gene_status)

#relevel the genes as normal as reference factor
applyrevel<-function(combined){
  
  col_names<-colnames(combined)[7:ncol(combined)]
  for(i in col_names){
    #combined$i<-as.factor(combined$i)
    combined[,i]<-as.factor(combined[,i])
    combined[,i]<-relevel(combined[,i],ref="Normal")
    #combined$i<-relevel(fctocombined$i,ref="Normal")
  }
  combined
}
combined<-applyrevel(combined)
######################################################
################ Univariate analysis ################
#####################################################

library(survival)
kmsurvo<-Surv(combined$Rectime,combined$Censor_status)
applycox<-function(combined){
  models<-list()
  col_names<-colnames(combined)[7:ncol(combined)]
  for(i in col_names){
    
    fit<-coxph(kmsurvo~factor(combined[,i]),data=combined)
    tss<-summary(fit)
    coefs<-c(tss$coefficients[c(1,3,5)])
    models[[i]]=coefs
  }
  
  final_mode<-as.data.frame(models) 
  final_model=t(final_mode)
  colnames(final_model)<-c("coef","exp.coef","p")
  
  as.data.frame(final_model)
}
fs<-applycox(combined)

#class(fs)
fs<-fs%>%mutate(gene=rownames(.))
fs
######################################################
################ Multivariate Analysis################
#####################################################

fitt<-coxph(kmsurvo~.,data=combined[,7:ncol(combined)])
fitt


######################################################
############### Survival Curve Plotting ##############
#####################################################

library(survival)
library(survminer)


for(i in 7:ncol(combined))
{
  kmsurvo <- Surv(combined$Rectime, combined$Censor_status)
  sfit <- survfit(Surv(combined$Rectime, combined$Censor_status) ~ (combined[,i]), data = combined)
  cox_fit <- coxph(Surv(combined$Rectime, combined$Censor_status) ~ (combined[,i]), data = combined)
  
  HR <- exp(coef(cox_fit))[[1]]
  
  ggsurv <- ggsurvplot(
    sfit,
    size = 1.5,
    #pval = TRUE,
    #pval.method = TRUE,
    #pval.coord = c(0, 0.05),
    legend.title = "",
    legend.labs = c("Normal", "Altered"),
    xlab = "Time in days",
    font.x = c(18, "bold.italic", "darkred"),
    font.y = c(18, "bold.italic", "darkred"),
    legend = c(0.75, 0.75),
    title = colnames(combined)[i],
    main = "Survival curve with Hazard Ratio",
    font.main = c(22, "bold", "darkblue"),
    font.tickslab = c(16, "plain", "darkgreen")
  )
  
  ggsurv$plot <- ggsurv$plot +
    theme(legend.text = element_text(size = 18, color = "black", face = "bold")) +
    annotate("text", x = 300, y = 0.25, 
             label = sprintf("HR = %.2f",HR), size = 6, color = "darkblue", fontface = "bold")
  
  print(ggsurv)
  
  ggsave(
    filename = paste0(colnames(combined)[i], ".png"),
    plot = ggsurv$plot,
    width = 7,
    height = 6,
    dpi = 300
  )
}
