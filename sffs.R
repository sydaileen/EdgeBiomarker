# Sequential Forward Floating Selection algorithm
# Usage: sffs_result=sffs(file,group,maxIter,kfold)

# file: l by n matrix of edge data, where l is the number
#		of edge features, n is the number of samples;
#		rownames are edge feature names taking the form:
#		"<name of node1>_<name of node2>_<class i>", 
#		where name of node1 is prior to name of node2 
#		in lexicographic order
# group: n by 1 matrix labeling the sample classes: 1,2,3,...,K;
#		rownames are sample names;
# sffs_result:a list containing two items
#			ft_seq: a list containing indexes of each feature set
#			J_seq:	a vector of cross validation accuracy for 
#					each selected feature set
# Created by Yidi Sun (sunyidi@picb.ac.cn)

load_packages=function(){
	a = function(){
		if ( require (e1071))
			return()
		else 
			install.packages("e1071")
		}
		
	b = function(){
		if ( require (ggplot2))
			return()
		else
			install.packages("ggplot2")
		}
	a()
	b()
}
load_packages()


sffs=function(file,group,maxIter=100,kfold=10){

	data=read.table(file,header=TRUE,row.names=1)
	labl=read.table(group,header=TRUE,row.names=1)
	# parameter adjustment
	if(kfold>dim(data)[2] || kfold==0){
		kfold=dim(data)[2]
		nrep=1
	}

	# data normalization for each feature (optional)
	data=t(scale(t(data)))

	iter=0
	universe=1:dim(data)[1]		# all features 
	ft_seq=list()						# recording optimal feature sets
	J_seq=rep(0,length=dim(data)[1])	# optimal scores of ft_seq
	current=c()		# current feature set to be considered

	require(e1071)
	# main body
	while(iter<maxIter){
		iter=iter+1
		# forwarding algorithm
		print(paste('iter=',iter,',forwarding...\n',sep=""),quote=FALSE)
		pool=universe[!(universe %in% current)]
		if(length(pool)==0){
			break
		}
		J_pool=c()
		for (i_p in 1:length(pool)){
			temp=svmCV(data[c(current,pool[i_p]),],labl,kfold)
			J_pool[i_p]=mean(temp)
		}
		J_max=max(J_pool)
		idx=which(J_pool==max(J_pool))
		if(J_seq[length(current)+1]<J_max){
			ft_seq[[length(current)+1]]=c(current,pool[idx])
			J_seq[length(current)+1]=J_max
			current=c(current,pool[idx])
		}else{
			current=ft_seq[[length(current)+1]]
		}
		print(paste('features=',length(current),";",'J_seq=',
			J_seq[length(current)],sep=""),quote=FALSE)

		# floating algorithm
		print(paste('iter=',iter,',floating...\n',sep=""),quote=FALSE)
		while(length(current)>2){
			J_ft=rep(0,length(current))
			for (i_f in 1:length(current)){
				ft_tmp=current[-i_f]
				temp=svmCV(data[ft_tmp,],labl,kfold)
				J_ft[i_f]=mean(temp)
			}
			J_max=max(J_ft)
			idx=which(J_ft==max(J_ft))
			if(J_seq[length(current)-1]<J_max){
				J_seq[length(current)-1]=J_max
				ft_seq[[length(current)-1]]=current[-idx]
				current=current[-idx]
				print(paste('features=',length(current),";",'J_seq=',
					J_seq[length(current)],sep=""),quote=FALSE)
			}else{
				break
			}
		}
	}
	J_seq=J_seq[J_seq!=0]
	# line plot
	colclass<-rep("numeric",2)
	colnames<-c("features","accuracy")
	acc<-read.table(text="",colClasses = colclass,col.names=colnames)
	for (k in 1:length(ft_seq)){
		acc[k,1]=length(ft_seq[[1]][[k]])
		acc[k,2]=J_seq[k]
	}
	require(ggplot2)
	pdf("sffs_result.pdf")
	ggplot(data=acc, aes(x=features, y=accuracy) + geom_line() + geom_point(size=1) 
		+ xlab("Number of Features") + ylab("Accuracy") 
		+ ggtitle("Cross validation accuracy")
		+ theme_bw()+theme(plot.title = element_text(hjust = 0.45))
	dev.off()
	sffs_result=list(ft_seq,J_seq)
	sffs_result
}

svmCV=function(data,labl,kfold=10){
	db=cbind(t(data),labl)
	colnames(db)[dim(db)[2]]="grp"
	sv=svm(grp~.,data=db,cross=kfold,type='C-classification',kernel='sigmoid') 
	acc=mean(sv$accuracies)
	acc
}