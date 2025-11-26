library(DESeq2)

### EDIT THE FOLLOWING ######
main_folder<-"MCF7" 
file_counts<-"selected_data.tsv"
file_metadata<-"selected_metadata.tsv"


####################

timestamp<-sub("selected_data_","",file_counts)


counts <- read.delim(file.path(main_folder,file_counts),sep="\t",row.names = 1)
metadata <- read.delim(file.path(main_folderfile_metadata),sep="\t")

counts <- counts[rowSums(counts) > 10, ] #remove low count genes before normalization

dds <- DESeqDataSetFromMatrix(
	countData = counts,
	colData = metadata,
	design = ~GROUP
)
dds <- DESeq(
	dds
)
res <- results(dds, contrast=c("GROUP","G1","G2"))
#res <- lfcShrink(dds, coef=paste(sep="_","GROUP","G2","vs","G1"), type="apeglm")

res_filtered<-res[complete.cases(res), ]
write.table(res_filtered, file = file.path(main_folder,paste0("de_data",timestamp)), sep = "\t")

res_filtered <- res_filtered[res_filtered$padj < 0.05 & abs(res_filtered$log2FoldChange) > 1.5 ,]
write.table(res_filtered, file = file.path(main_folder,paste0("de_data_subset",timestamp)), sep = "\t")

