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

groupCapture = function(input,g=1,p=5,lwd=2,pch=1,out=c("data","samples","bool"),data=NULL) {
	res=list()
	plot(input[,1],input[,2],pch=pch)
	for (i in 1:g) {
		click <- locator(p)
		clickX=c(click$x,click$x[1])
		clickY=c(click$y,click$y[1])
		lines(clickX,clickY,col=i,lwd=lwd)
		d=data.frame(clickX,clickY)
		sel=in.out(as.matrix(d),as.matrix(input))
		if (out=="data") {
			if(exists("data")) { 
				res[[i]]=data[,sel]
			} else {
				cat("no data provided")
				break()
			}				
		}
		if(out=="samples") {
			res[[i]]=rownames(input)[sel]
		}
		if(out=="bool") {
			res[[i]]=sel
		}
	}
	res
}
