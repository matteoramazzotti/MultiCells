library(httr)
urla="https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
annot <- data.table::fread(urla, verbose=F,header=T, quote="", stringsAsFactors=F, data.table=F)
out=list()
#this function downloads counts, annotations and sample names in the form of a list of data frames
RNAseq_autodown=function(id) {
	urlf="https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&acc=XXX&file=XXX_raw_counts_GRCh38.p13_NCBI.tsv.gz"
	path=gsub("XXX",id,urlf)
	if(isTRUE(!sapply(path, http_error))) {
		tbl=as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
		urln="https://ftp.ncbi.nlm.nih.gov/geo/series/YYY/XXX/soft/XXX_family.soft.gz"
		urln=gsub("XXX",id,urln)
		urln=gsub("YYY",gsub("\\d\\d\\d$","nnn",id,perl=T),urln)
		urln=paste0("wget -q -O- ",urln," | zcat | grep Sample_title | perl -ne 's/!Sample_title = //;print'")
		names=read.delim(pipe(urln),header=F)[,1]
		out=list()
		out$counts=as.matrix(tbl)
		out$names=names
		out
	}
	else {
		cat("Broken link\n")
	}
}
