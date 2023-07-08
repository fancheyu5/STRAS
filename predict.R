args<-commandArgs(trailingOnly = TRUE)

library(stringr)
library(dplyr)
library(tidyr)
library(caret)


str2max_num<-function(string){
	number<-as.numeric(unlist(str_split(string,";")))
	return(max(number))
}

rffit <- readRDS('rffit.2.1.3.rds')

dataset<-read.table(file=args[1],sep='\t',header=F,quote="")
colnames(dataset)<-c("chr","start","end","patient_motif","patient_copynumber","genetype",
	                  "genename","protien_coding","exon","UTR","ref_motif",
	                  "ref_period","mean_1000g_50Han","sd_1000g_50Han","max_1000g_50Han",
	                  "mean_1000g_100","sd_1000g_100","max_1000g_100","TADboundry",
	                  "CTCF","TFBS","gene_pLI","HPO")
dataset<-mutate(dataset,ID = str_c(chr,start,end,patient_copynumber,sep="_"))
data_for_predict<-mutate(dataset,
	                 num_mean_1000g_100=sapply(dataset$mean_1000g_100,str2max_num),
	                 num_sd_1000g_100=sapply(dataset$sd_1000g_100,str2max_num),
	                 num_max_1000g_100=sapply(dataset$max_1000g_100,str2max_num),
	                 num_mean_1000g_50Han=sapply(dataset$mean_1000g_50Han,str2max_num),
	                 num_sd_1000g_50Han=sapply(dataset$sd_1000g_50Han,str2max_num),
	                 num_max_1000g_50Han=sapply(dataset$max_1000g_50Han,str2max_num),
                         expansion = ifelse(patient_copynumber > num_mean_1000g_100,T,F),
	                 expansion_rate = (patient_copynumber-num_max_1000g_100)/num_max_1000g_100, 
	                 expansion_length =  patient_copynumber-num_mean_1000g_100, 
	                 peroid_length = str_length(patient_motif),
	                 sd = num_sd_1000g_100,                             
	                 genetic = ifelse(genetype != "intergenic",T,F),
	                 protien_coding =ifelse(protien_coding != "-",T,F),
	                 exon = ifelse(exon != "-",T,F),
	                 CTCF = ifelse(CTCF != "-",T,F),
	                 TFBS = ifelse(TFBS != "-",T,F),
	                 TADboundry = ifelse(TADboundry != "-",T,F),
	                 popdiff = num_sd_1000g_100-num_sd_1000g_50Han,
	                 known_disease_associated_gene = ifelse(HPO != "-",T,F))
predicted<-na.omit(subset(data_for_predict,peroid_length>2))
predictY<-predict(rffit, newdata = predicted, type = "prob")

score_bed<-predicted[,c("ID","expansion_rate")]
score_bed$score<-predictY[,1]
colnames(score_bed)<-c("ID","expansion_rate","score")

scored_data<-left_join(dataset,score_bed,by="ID")
sorted_data<-scored_data[order(scored_data$score,na.last=T,decreasing=T),]
write.table(sorted_data,file=args[2],row.names=F,col.names=T,sep='\t')
