##### GIVEN A SELECTION OF SAMPLES FROM seq_exp_log_XXX_ok
coll=read.delim("seq_exp_log_MCF7_ok",header=F)
#coll[,1]=unique(sub("^ +","",perl=T,coll[,1]))

#### DOWNLOAD AND ASSEMBLE THE SELECTED SEQ DATA
#### CAVEAT: in case samples are in two different GPL, it probably fails. 
#### In this case, remove the GSE (no other solution found, 3 days to get here)

library(httr)
#generalized url of the count data.
urlf="https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&acc=XXX&file=XXX_raw_counts_GRCh38.p13_NCBI.tsv.gz"
#generalized url of the sample annotations
urln="https://ftp.ncbi.nlm.nih.gov/geo/series/YYY/XXX/soft/XXX_family.soft.gz"

#coll[,1]=unique(sub(" GSE","GSE",coll[,1]))
#makes GSE unique
GSEs=unique(gsub(" ","",coll[,1]))
data=list()
data_colnames="Gene"

#for each GSE
#in case some GSE breaks the download, jump to the next gse by chencing the start below
start=1
for (i in start:length(GSEs)) {
	#catch selected sample names in GSE 
	samples=coll[grep(GSEs[i],coll[,1]),2]
	cat(paste0("--- ",i,"/",length(GSEs),") ", GSEs[i])," Sample/s: ",paste0(samples,sep="|"),"---\n",sep="")
	#adjust the url for downloading count data
	path=gsub("XXX",GSEs[i],urlf)
	#download counts
	tbl=as.data.frame(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
	#adjust the url for downloading actual names of samples (will be matched with the samples above)
	url=gsub("XXX",GSEs[i],urln)
	url=gsub("YYY",paste0(substr(GSEs[i],1,nchar(GSEs[i])-3),"nnn"),url)
	#prepare the link to download sample names
	url=paste0("wget -q -O- ",url," | zcat | grep Sample_title | perl -ne 's/!Sample_title = //;print'")
	cat(paste0("URL: ",url,"\n"))
	#download sample names
	names=read.delim(pipe(url),header=F)[,1]
	if(i==1) {
		#adds gene names as the first column to initialize the data data.frame
		col_one=tbl[,1]
	}
	#samples are the selected ones, they should match those obtained form NCBI
	found=match(samples,names)
	       #if all saple are matched                 AND  this is needed for malformed tbl: the whole exp must be removed... 
	if(sum(found %in% 1:dim(tbl)[2]) == length(found) && length(tbl)>length(found)) {
		#this +1 is because the tbl matrix has gene names as first row, so the matched columns are to be shifted by 1
		TMP=tbl[,found+1]
	} else {
		cat("Sample/s not found\n")
		next
	}
	#print(data.frame(samples,found,names[found],colnames(TMP)))
	data[[i]]=data.frame(TMP)
	#store sample names that will be used as column names for the main data table (whitespaces turns into "_")
	data_colnames=gsub(" ","_",paste0(GSEs[i],"_",names[match(samples,names)]))
	if(class(data[[i]])=="integer") {
		#this happens when only 1 sample is selected fomr the downloaded data (so it is a vector, not a data.frame) 
		names(data[[i]])=data_colnames
		cat("Current table column is 1\n")
	} else {
		#more than one sample, so this is a data.frame
		colnames(data[[i]])=data_colnames
		cat("Current table column is",dim(data[[i]])[2],"\n")
	}
	print(head(data[[i]]))
}

#this creates the data matrix for this experiment
#seq_matrix = data[[1]]
seq_matrix = data.frame(col_one)
for (i in 1:length(data)) {
	cat("dataset",i)
	if (!is.null(dim(data[[i]]))) {
		if (is.data.frame(data[[i]])) {
			seq_matrix=cbind(seq_matrix,data[[i]])
			cat(" is a data frame of length", dim(data[[i]])[2],"\n")
		} else {
			seq_matrix=cbind(seq_matrix,as.numeric(data[[i]]))
			cat(" is numeric\n")
		}
	} else {
		cat(" is missing\n")
	}
}
rownames(seq_matrix)=col_one
seq_matrix=seq_matrix[,-1]

save(file="MCF7.RData",GSEs,seq_matrix)
load(file="MCF7.RData")


#### download annotations

urla="https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
annot <- data.table::fread(urla, verbose=F,header=T, quote="", stringsAsFactors=F, data.table=F)[,1:3]
#or
#system("wget https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz")
#annot=read.delim("../../Human.GRCh38.p13.annot.tsv.gz")[,1:3]
