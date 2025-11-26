library(DESeq2)
library(getopt)

spec = matrix(c(
  'verbose', 'v', 2, "integer", # logical, integer, double, complex, character.
  'help'   , 'h', 0, "logical",
	'counts_file','c',1,"character",
	'metadata_file', 'm',1,"character",
	'main_folder', 'M',1,"character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

main_folder<-opt$main_folder
file_counts<-opt$counts_file
file_metadata<-opt$metadata_file


timestamp<-sub(paste0(main_folder,"/","selected_data_"),"",file_counts)

counts <- read.delim(file_counts,sep="\t",row.names = 1)
metadata <- read.delim(file_metadata,sep="\t")

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
write.table(res_filtered, file = file.path(main_folder,paste0("de_data_",timestamp)), sep = "\t")

res_filtered <- res_filtered[res_filtered$padj < 0.05 & abs(res_filtered$log2FoldChange) > 1.5 ,]
write.table(res_filtered, file = file.path(main_folder,paste0("de_data_subset_",timestamp)), sep = "\t")

