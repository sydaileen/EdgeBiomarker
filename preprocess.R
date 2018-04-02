# Preprocess data before network construction
# Usage: node_data = preprocess(path, data, outfile)
# data: m by n matrix of node data, where m is the number
#       of node features, n is the number of samples; 
#		rownames are feature names; colnames are sample names

preprocess = function (path, data, method = "zscore", outfile="Proprocess.xls"){
	setwd(path)
	db=read.table(data,header=TRUE,row.names=1)
	# optional 
	db=log2(db+1)

	pdf("Data_raw.pdf")
	boxplot(db)
	dev.off()

	# z-score scale
	if (method == "zscore"){
     	db=t(scale(t(db)))  
  	}

  	# median scale
	if (method == "median"){
     	db=apply(db,1,function(x){x/median(x, na.rm=T)})  
  	}

  	# mean scale
	if (method == "mean"){
     	db=apply(db,1,function(x){x/mean(x, na.rm=T)})  
  	}
	
	pdf("Data_norm.pdf")
	boxplot(db)
	dev.off()
	
	write.table(db, file=outfile, quote=FALSE,sep="\t")
	db
}
