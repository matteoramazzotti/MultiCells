library(mgcv)
#example data
#data=data.frame(matrix(rnorm(10000),nrow=100),row.names = paste0("G", 1:100))
#colnames(data)=paste0("S", 1:100)

#how to perform dimensionality reduction to create
#orthogonal coordinated for plotting
#using R-base PCA
#test=prcomp(data)$rotation[,1:2]
#using Rtsne
#library(Rtsne)
#test=Rtsne(t(as.matrix(data)), dims = 2, perplexity=4, verbose=TRUE, max_iter = 500,check_duplicates = FALSE)
#using Nonmetric Multidimensional Scaling (from the vegan package) 
#test=metaMDS(as.matrix(data))$points

#input = 2D matrix with orthogona coordinates
#g = number of groups to form
#p = number of vertex of the selection polygon
#out = "data" emits a list of matrices (subset of the main one) containing data of the selcted samples (the data= input must be present)
#out = "samples" emits a list of sample names per group
#out = "bool" emits a list of true/false to be used in external selection

groupCapture = function(
	input,
	data=NULL,
	xlab=c(""),
	ylab=c(""),
	g=1,
	p=5,
	lwd=2,
	pch=1,
	out=c("data","samples","bool","df"),
	auto_scale=TRUE,
	range=c(-100,100),
	labels_on = TRUE,
	labels = c("")
) {
	res=list()

	if (auto_scale) {
		x_min <- min(input[,1])
		x_max <- max(input[,1])
		y_min <- min(input[,2])
		y_max <- max(input[,2])
		max <- abs(round(max(x_max,y_max) + 10))
		min <- abs(round(min(x_min,y_min) - 10))
		range=c(-min,max)
	}
	plot(
		input[,1], input[,2],
		xlim = range,
		ylim = range,
		xlab = xlab,
		ylab = ylab,
		pch = pch,
		col = "steelblue"
	)
	text(
		pca$x[,1], 
		pca$x[,2], 
		labels = labels,
		pos = 3, 
		cex = 0.3
	)
	for (i in 1:g) {
		click <- locator(p,type="l")
		clickX=c(click$x,click$x[1])
		clickY=c(click$y,click$y[1])
		lines(clickX,clickY,col=i,lwd=lwd)
		d=data.frame(clickX,clickY)
		sel=in.out(as.matrix(d),as.matrix(input))
		if (out[1]=="data" || out[1]=="df") {
			if(exists("data")) {
				ind<-colnames(data)[sel]
				res[[i]]=data.frame(data[,ind])
			} else {
				cat("no data provided")
				break()
			}				
		}
		if(out[1]=="samples") {
			res[[i]]=rownames(input)[sel]
		}
		if(out[1]=="bool") {
			res[[i]]=sel
		}

	}
	if (out[1] == "df") {
		tmp = res[[1]]
		tmp2 = data.frame(SAMPLE=colnames(tmp),GROUP=sapply(1:length(colnames(tmp)),function(x){paste0("G","1")}))

		if (g>1) {
			for (i in 2:g) {
				tmp = cbind(tmp,res[[i]])
				tmp2 = rbind(tmp2,data.frame(SAMPLE=colnames(res[[i]]),GROUP=sapply(1:length(colnames(res[[i]])),function(x){paste0("G",i)})))
			}
		}
		res = list(tmp,tmp2)
	}
	return(res)
}
