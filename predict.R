args<-commandArgs(trailingOnly = TRUE)

library(stringr)
library(dplyr)
library(tidyr)
library(caret)


str2max_num<-function(string){
	number<-as.numeric(unlist(str_split(string,";")))
	return(max(number))
}

rffit <- readRDS('rffit.1.1.rds')

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
                     expansion = if_else(patient_copynumber > num_mean_1000g_100,T,F),
	                 expansion_rate = (patient_copynumber-num_max_1000g_100)/num_max_1000g_100, 
	                 expansion_length =  patient_copynumber-num_mean_1000g_100, 
	                 peroid_length = str_length(patient_motif),
	                 sd = num_sd_1000g_100,                             
	                 genetic = if_else(genetype != "intergenic",T,F),
	                 protien_coding =if_else(protien_coding != "-",T,F),
	                 exon = if_else(exon != "-",T,F),
	                 CTCF = if_else(CTCF != "-",T,F),
	                 TFBS = if_else(TFBS != "-",T,F),
	                 TADboundry = if_else(TADboundry != "-",T,F),
	                 popdiff = num_sd_1000g_100-num_sd_1000g_50Han,
	                 known_disease_associated_gene = if_else(HPO != "-",T,F))
predicted<-na.omit(subset(data_for_predict,peroid_length>2))
predictY<-predict(rffit, newdata = predicted, type = "prob")

score_bed<-as.data.frame(predicted$ID)
score_bed$score<-predictY[,1]
colnames(score_bed)<-c("ID","score")

scored_data<-left_join(dataset,score_bed,by="ID")
sorted_data<-scored_data[order(scored_data$score,na.last=T,decreasing=T),]
write.table(sorted_data,file=args[2],row.names=F,col.names=T,sep='\t')
