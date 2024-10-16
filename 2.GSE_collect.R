##### FUNCTIONS TO DOWNLOAD SAMPLE INDEXES
#this procedure saves results in the file "seq_exp_log"

#delete existing log files
system("unlink seq_exp_log")

source("RNAseq_autodown.R")
seq_exp=read.delim("GSE_seq")[,1]
seq_dataset=NULL
seq_track=NULL
for (i in 1:length(seq_exp)) {
	if(!seq_exp[i] %in% seq_track){
		cat("Getting",i,"/",length(seq_exp),seq_exp[i],"\n")
		seq_dataset=RNAseq_autodown(seq_exp[i])
		cat(paste0(seq_exp[[i]],"\t",seq_dataset$name,"\n"),file="seq_exp_log",append=T)
		seq_track=c(seq_track,seq_exp[i])
	} else {
		cat(seq_exp[i],"already downloaded\n")
	}
}

#the seq_exp_log file is filtered to select rows (experiments) matching the cell name MCF7 
system("sort seq_exp_log | uniq | grep -i MCF7 | perl -ne 's/^ //;print' > seq_exp_log_MCF7")

#this strips white specaes at the beginning of the line
system("perl -ne 's/^ //;print' seq_exp_log_MCF7 > seq_exp_log_MCF7")

#manual selection of experiments based on the "seq_exp_log_MCF7" file
#the aim is to create a two column data frame with
#GSExxxx    sample_name
#named seq_exp_log_MCF7_ok
