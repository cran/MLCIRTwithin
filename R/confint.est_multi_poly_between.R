confint.est_multi_poly_between <-function(object, parm, level=0.95, ...){
	
# preliminaries
	parm = NULL
	out = object
# print output
	if(is.null(out$seTh)) stop("Use option out_se=TRUE in est_multi_poly_within()")
	za = qnorm(1-(1-level)/2)	
    cat("\nConfidence interval for abilities:\n")
    L1Th = out$Th-za*out$seTh
    L2Th = out$Th+za*out$seTh
    Tmp = cbind(L1Th,L2Th)
    ind = NULL
    for(j in 1:ncol(out$Th)) ind = c(ind,j,j+ncol(out$Th))
	Tmp = matrix(Tmp[,ind],nrow(out$Th),2*ncol(out$Th))
    dimnames(Tmp) = list(dimension=1:nrow(out$Th),class=(1:ncol(out$Th))%x%c(1,1))
    print(round(Tmp,4))    
    cat("\nConfidence interval for the item parameters:\n")
    Tmp = out$Bec-za*out$seBec
    for(j in 1:ncol(out$Bec)) colnames(Tmp)[j] = paste("beta",j,sep="")
    L1Items = cbind(gamma=out$ga1c-za*out$segac,Tmp)
    Tmp = out$Bec+za*out$seBec
    for(j in 1:ncol(out$Bec)) colnames(Tmp)[j] = paste("beta",j,sep="")
    L2Items = cbind(gamma=out$gac+za*out$segac,Tmp)
    Tmp = cbind(L1Items,L2Items)
    ind = NULL
    for(j in 1:ncol(L1Items)) ind = c(ind,j,j+ncol(L1Items))
    Tmp = Tmp[,ind]
    print(round(Tmp,4))              
    cat("\nConfidence interval for the regression coefficients:\n")
    L1De = out$De-za*out$seDe
    L2De = out$De+za*out$seDe
   	Tmp = cbind(L1De,L2De)
    ind = NULL
    for(j in 1:ncol(out$De)) ind = c(ind,j,j+ncol(out$De))
    Tmp = Tmp[,ind]
    names = dimnames(Tmp)
    dimnames(Tmp) = list(names[[1]],logit=names[[2]])
    print(round(Tmp,4))
    cat("\n")
# output
	out = list(L1Th=L1Th,L2Th=L2Th,L1Items=L1Items,L2Items=L2Items,L1De=L1De,L2De=L2De)    
	
}