library(clusterProfiler)
# library(enrichplot)
library(ggplot2)
# library(msigdbr)
library(getopt)

spec = matrix(c(
  'verbose', 'v', 2, "integer", # logical, integer, double, complex, character.
  'help'   , 'h', 0, "logical",
  'genes'  , 'g', 1, "character",
  'gmt'   , 'G', 1, "character",
  'filter' , 'f', 2, "character",
	'out-prefix', 'o', 1, "character",
	'out-dir', 'O', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


dataTable <- read.delim(opt$genes,sep="\t")
myInput <- as.vector(rownames(dataTable))

custom_gmt <- opt$gmt
gene_sets <- clusterProfiler::read.gmt(custom_gmt)
# gene_sets$term <- gsub("^GOBP_", "", gene_sets$term)
# gene_sets$term <- gsub("_", " ", gene_sets$term)           
# gene_sets$term <- tools::toTitleCase(tolower(gene_sets$term))


# m_df <- msigdbr(species = "Homo sapiens")
# m_t2g <- msigdbr(species = "Homo sapiens", collection = "C5") %>% 
#   dplyr::select(gs_name, gene_symbol)


enriched <- clusterProfiler::enricher(
	gene = myInput,
	TERM2GENE = gene_sets,
)
myDF <- data.frame(enriched$ID,enriched$Count,enriched$GeneRatio,enriched$pvalue,enriched$p.adjust)
write.table(myDF,paste(opt['out-dir'],"/","cp_results_",opt['out-prefix'],".tsv",sep=""), sep="\t",row.names=FALSE,quote=FALSE)

df_top<-myDF[order(myDF$enriched.p.adjust,decreasing = FALSE),]
df_top<-head(df_top,20)

df_top$enriched.ID <- factor(
	df_top$enriched.ID,
  levels = rev(df_top$enriched.ID)
)


ggplot(
	df_top,
	aes(
		x = enriched.Count,
    y = enriched.ID,
    fill = enriched.p.adjust)
) +
geom_col() +
scale_fill_gradientn(
	# high = "red",
	# low = "blue",
	colours = c( "#0000FF", "#808080", "#FF8800"),
	name = "adj. p-value",
	trans = "reverse"
) +
theme_minimal(base_size = 14) +
theme(
	plot.background = element_rect(fill="white"),
	axis.text.y = element_text(size = 5,face="bold"),  # shrink labels if many categories
	axis.text.x = element_text(size = 5,face="bold"),  
	axis.title.y = element_blank(),
	axis.title.x = element_text(size = 8),
	plot.margin = margin(10, 10, 10, 10),   # adjust margins
	legend.text = element_text( size=5),
	legend.title = element_text(size=8,hjust=1),
	legend.margin = margin(10, 5, 10, 10)
) +
labs(
	x = "Gene Count",
	title = ""
)
ggsave(paste(opt['out-dir'],"/plots/","barplot_top_",opt['out-prefix'],".png",sep=""),
	width = 1920,
  height = 1080,
  units = "px",
  dpi = 300,)