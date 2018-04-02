# Prognostic feature selection by LASSO cox model 
# Usage: lassotable = lassocox (file, clinical, nfold, outfile)
# file: p by n matrix of data, where p is the number
#       of features, n is the number of samples; 
#       rownames are feature names; colnames are sample names
# clinical: n by q matrix of clinical data, where n is 
#           the number of samples, q is the clinical informtion,
#           including at least survival time and survival event
# lasso: r by 2 lasso coeficients matrix, r is the number of 
#         selected prognostic feature values,the first column is 
#         the feature names, the second column is the 
#         coefficient of each feature 
# Created by Yidi Sun (sunyidi@picb.ac.cn)

# load required packages
load_packages=function(){
  a = function(){
    if ( require (glmnet))
      return()
    else 
      install.packages("glmnet")
    }
    
  b = function(){
    if ( require (survival))
      return()
    else
      install.packages("survival")
    }
    
  a()
  b()
}
load_packages()

# LASSO regression
require(glmnet)
require(survival)

lassocox = function(file, clinical, nfold=5, outfile="prog_features.xls"){
  data=read.table(file, header=TRUE, row.names=1)
  clin=read.table(clinical, header=TRUE, row.names=1)
  y = Surv(clin$surv_year,clin$vital_status)
  feature.p = c()
  # transpose the data set 
  data=t(data)
  for (i in 1:ncol(data)){
    cox = coxph(y ~ data[,i],method="breslow")
    feature.p[i] = summary(cox)$logtest["pvalue"]
    names(feature.p)[i] = colnames(data)[i]
  }
  feature.p.sig = feature.p[feature.p<=0.05]
  dbnew = data[,colnames(data) %in% names(feature.p.sig)]

  fit.cv  =  cv.glmnet(x=dbnew,y=y,family="cox",alpha=1,standardize=FALSE,nfolds=5) 
  lambda  =  fit.cv$lambda.min
  final  =  fit.cv$glmnet.fit
  coef.fit  =  coef(final,s=lambda)[2:(ncol(dbnew)+1)]
  cols.include  =  which(abs(coef.fit) > 0)
  feature.name = colnames(dbnew)[cols.include]
  lasso = data.frame(cbind(feature.name,coef.fit[cols.include]))
  write.table(lasso,file=outfile,sep="\t",quote=FALSE)
  lasso
}
