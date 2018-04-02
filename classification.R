# demonstrate the classification effect of selected features
# Usage: features = classification (file, group, sffs_result, 
# method="both",outfile="features.txt")
# data: p by n matrix of data, where p is the number
#       of features, n is the number of samples; 
#		rownames are feature names; colnames are sample names
# sffs_result: a list containing two items
#			ft_seq: a list containing indexes of each feature set
#			J_seq:	a vector of cross validation accuracy for 
#					each selected feature set
# Created by Yidi Sun (sunyidi@picb.ac.cn)

# load required packages
load_packages=function(){
	a = function(){
		if ( require (ggplot2))
			return()
		else 
			install.packages("ggplot2")
		}
		
	b = function(){
		if ( require (gplots))
			return()
		else
			install.packages("gplots")
		}
	
	c = function(){
		if (require (pheatmap))
			return()
		else 
			install.packages("pheatmap")
		}
		
	a()
	b()
	c()
}
load_packages()

classification = function (file, group, sffs_result, method="both",outfile="features.txt"){
	data=read.table(file, header=T,row.names=1)
	labl=read.table(group,header=TRUE,row.names=1)

	ft_seq=sffs_result[[1]]
	J_seq=sffs_result[[2]]
	tmpid=which(J_seq == max(J_seq))
	idx=ft_seq[[tmpid]]
	featable=data[idx,]
	write.table(featable, file=outfile,sep="\t",quote=FALSE)

	# pca plot
	plot1=pca(featable, group=labl[,1])

	#hca plot
	plot2=heatmap(featable, group=labl)

	if(method=="pca"){
		pdf("PCA.pdf")
		plot1
		dev.off()
	}
	if(method=="hca"){
		pdf("HCA.pdf")
		plot2
		dev.off
	}
	if(method=="both"){
		pdf("PCA.pdf")
		plot1
		dev.off()
		pdf("HCA.pdf")
		plot2
		dev.off()
	}
}

pca=function(data,group=1,title = NULL,legend="right",text=FALSE,xlim=NULL,ylim=NULL,size.point=3,size.text=2){
	group=as.factor(group)
	if (length(group)==1) legend.position = "none" else legend.position = legend
	#calculation
	pr=prcomp(data,cor=T)  
	Comp=pr$sdev
	Comp.var=Comp*Comp

	Comp.1=paste("Comp.1(",round(Comp.var[1]/sum(Comp.var)*100,2),"%)",sep="")
	Comp.2=paste("Comp.2(",round(Comp.var[2]/sum(Comp.var)*100,2),"%)",sep="")
	pcomp=pr$x
	
	#plot
    label=as.character(rownames(pcomp))
	require(ggplot2)
	p=ggplot(as.data.frame(pcomp),aes(PC1,PC2))
	coord = coord_fixed(ratio =Comp.var[2]/Comp.var[1], xlim = xlim, ylim = ylim)
	theme1=theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
	geom = geom_point(size=size.point,shape=16,aes(colour=group))
	theme2=theme(legend.position=legend.position)+theme(plot.title = element_text(size = 20))+theme(axis.title = element_text(size = 12))
	lab=labs(x=Comp.1,y=Comp.2,title=title,size=10)
	if (text==FALSE) plot=p+coord+theme1+geom+theme2+lab
	else  {point.text= geom_text(label=label,size=size.text,hjust=1,vjust=1);
			plot=p+coord+theme1+geom+theme2+lab+point.text}
	    
	plot		
}

heatmap=function(data,scale="row",k_rows=1,k_cols=1,file=NA,fontsize_row=2,fontsize_col=4,group=1,col = colorRampPalette(c("blue", "white", "red"))(297),main=""){
    mat=data
    require(pheatmap)
	require(gplots)   #for greenre(75))

	#heatmpa annotation
	if (length(group)==1) {
		annotation_col=NA
	} else {
		annotation_col=group
	}
   
   	pheatmap(mat,color=col,
         clustering_distance_row = "correlation",
         clustering_distance_col = "correlation",
         scale="row",cutree_rows=k_rows,cutree_cols=k_cols,
         fontsize=30,fontsize_row=8,fontsize.col=6,
         annotation_col=annotation_col,
         cell_width=8,cell_height=6,
         treeheight_row=200,treeheight_col=200,
         width=1000,height=6000,
         main=main)
}
