library(DESeq2)
source("groupCapture.R")

counts <- read.delim("matrix_symbol.tsv",sep="\t",row.names = 1)
# gsm_to_gse <- read.delim("gsm_to_gse_selected_no_singlets.tsv",sep="\t",header=FALSE )
gsm_to_gse <- read.delim("gsm_to_gse_selected.tsv",sep="\t",header=FALSE )
counts <- counts[rowSums(counts) > 10, ] #remove low count genes before normalization


sample_ids <- gsm_to_gse$V1
counts_subset <- counts[, colnames(counts) %in% sample_ids] #remove singlets


dds <- DESeqDataSetFromMatrix(
	countData = counts_subset,
	colData = DataFrame(row.names = colnames(counts_subset)),
	design = ~ 1
)

vst_data <- vst(dds, blind = TRUE) #vst normalization
vst_mat <- assay(vst_data)
v <- apply(vst_mat, 1, var)
vst_mat <- vst_mat[v > 0, ]

pca <- prcomp(t(vst_mat),scale=TRUE)
percentVar <- pca$sdev^2 / sum(pca$sdev^2)

# preparation for plot
match_indexes = match(colnames(vst_mat),gsm_to_gse[,1])
colnames_by_gse = gsm_to_gse[match_indexes,2] # labels for points


xlab = paste0("PC1: ", round(percentVar[1]*100, 1), "% variance") # label for x axis
ylab = paste0("PC2: ", round(percentVar[2]*100, 1), "% variance") # label for y axis


pca_df <- data.frame(
	PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

pca_df<-pca_df[complete.cases(pca_df), ]


# plot and select

res <- groupCapture(
	pca_df,
	data = counts_subset,
	g = 2,
	xlab = xlab,
	ylab = ylab,
	pch = 19,
	out = "df",
	auto_scale=FALSE,
	range = c(-200,200),
	labels = colnames_by_gse
)

metadata = data.frame(
	SAMPLE=res[[2]]$SAMPLE,
	SERIES=gsm_to_gse[match(res[[2]]$SAMPLE,gsm_to_gse[,1]),2],
	GROUP=res[[2]]$GROUP
)



# export 
write.table(res[[1]], file = "selected_data.tsv", sep = "\t")
write.table(metadata, file = "selected_metadata.tsv", sep = "\t",row.names=FALSE)