library(DESeq2)

counts <- read.delim("selected_data.tsv",sep="\t",row.names = 1)
metadata <- read.delim("selected_metadata.tsv",sep="\t")

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
write.table(res_filtered, file = "de_data.tsv", sep = "\t")

res_filtered <- res_filtered[res_filtered$padj < 0.05 & abs(res_filtered$log2FoldChange) > 1.5 ,]
write.table(res_filtered, file = "de_data_subset.tsv", sep = "\t")

