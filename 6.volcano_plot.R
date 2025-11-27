library(optparse)
library(ggplot2)

option_list <- list(
  make_option(
    c("-v", "--verbose"),
    type    = "integer",
    default = 0,
    help    = "verbose level (integer, optional)"
  ),

  make_option(
    c("-i", "--id"),
    type    = "character",
		default = NULL,
    help    = "analysis id (required)"
  ),
	make_option(
    c("-R", "--deseq2_results"),
    type = "character",
		default = NULL,
    help    = "tsv file holding DESeq2 results"
  ),
	make_option(
    c("-t", "--p_value_thr"),
    default = 0.05,
    help    = "Threshold for p-value (default: 0.05)"
  ),
	make_option(
    c("-T", "--log2_fc_thr"),
    default = 1.5,
    help    = "Threshold for log2FC (default: 1.5)"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

p_value_thr<-opt$p_value_thr
log2_fc_thr<-opt$log2_fc_thr



if (is.null(opt$deseq2_results) || is.null(opt$id)) {
	print("Missing -i or -R arguments")
	q()
} 

analysis_id<-c(opt$id)
files<-c(opt$deseq2_results)

for (file in files) {
	filename<-paste0(file,"_volcano")
	filename<-sub(".tsv","",filename)
	filename<-sub(analysis_id,"",filename)

	# pdf(file.path(".",paste(file,"volcano.pdf",sep="_")))
	df<-read.delim(file,sep="\t")
	df_clean<-df[complete.cases(df), ]
	df_clean<-df_clean[df_clean$padj > 0,]
	df_clean$diffexpressed <- "NO"
	df_clean$diffexpressed[df_clean$log2FoldChange > log2_fc_thr & df_clean$padj < p_value_thr] <- "UP"
	df_clean$diffexpressed[df_clean$log2FoldChange < -log2_fc_thr & df_clean$padj < p_value_thr] <- "DOWN"
	# df$delabel <- ifelse(df$X %in% head(df[order(df$pvalue), "X"], 30), df$X, NA)

	df_clean<-rbind(
		df_clean[df_clean$diffexpressed == "NO",],
		df_clean[df_clean$diffexpressed == "UP",],
		df_clean[df_clean$diffexpressed == "DOWN",]
	)

	color_values<-c()
	color_labels<-c()

	if (any(grepl("DOWN",df_clean$diffexpressed))) {
		color_values<-c(color_values,"#00AFBB")
		color_labels<-c(color_labels,"DOWN")
	}

	color_values<-c(color_values,"grey")
	color_labels<-c(color_labels,"NO")

	

	if (any(grepl("UP",df_clean$diffexpressed))) {
		color_values<-c(color_values,"#bb0c00")
		color_labels<-c(color_labels,"UP")
	}
	df_clean$neglogp <- -log10(df_clean$padj)
	CAP_Y <- 50                 # you can set 30–100 based on dataset
	df_clean$neglogp_cap <- pmin(df_clean$neglogp, CAP_Y)
	xlimit <- ceiling(max(c(
		abs(max(df_clean$log2FoldChange)),
		abs(min(df_clean$log2FoldChange))
	)) / 2.5) * 2.5

	ylimit <- -log10(min(df_clean$padj))
	# ylimit <- CAP_Y



	# pdf(file.path("analysis",analysis_id,"plots","volcano",filename))
	p<-ggplot2::ggplot(
		data = df_clean, 
		ggplot2::aes(
			x = .data[["log2FoldChange"]],
			y = -log10(.data[["padj"]]),
			# y = neglogp_cap,
			col = .data[["diffexpressed"]],
			# label = .data[["delabel"]]
		)
	) 	+
  	ggplot2::geom_vline(xintercept = c(-log2_fc_thr, log2_fc_thr), col = "gray", linetype = 'dashed') +
  	ggplot2::geom_hline(yintercept = -log10(p_value_thr), col = "gray", linetype = 'dashed') +
  	ggplot2::geom_point(size = 0.5) +
  	ggplot2::scale_color_manual(
			values = color_values, 
			labels = color_labels
		) +
  	ggplot2::coord_cartesian(ylim = c(0, ylimit), xlim = c(-xlimit, xlimit)) +
  	ggplot2::labs(
			color = 'Expression', #legend_title
      x = expression("log"[2]*"FC"),
			y = expression("-log"[10]*"p-value adj.")
		) +
 	 	# scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  	ggplot2::ggtitle('')  # Plot title
  	# geom_text_repel(max.overlaps = Inf) # To show all labels 
	# print(p)
	# dev.off()
	filepath<-file.path(analysis_id,"plots",filename)
	ggsave(
		filename = paste0(filepath,".png"),
		plot = p,
		dpi = 300,
		width = 1800, height = 2400, units = "px"
	)

	ggsave(
		filename = paste0(filepath,".tiff"),
		plot = p,
		dpi = 300,
		width = 1800, height = 2400, units = "px",
		compression = "lzw"
	)

	ggsave(
		filename = paste0(filepath,".pdf"),
		plot = p,
		width = 1800, height = 2400, units = "px"
	)

}



