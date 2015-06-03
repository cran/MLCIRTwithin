est_multi_poly_within <- function(S,yv=rep(1,ns),k1,k2,X=NULL,start=0,link=1,disc=0,difl=0,
                                  multi1=1:J,multi2=1:J,piv1=NULL,piv2=NULL,
                                  Phi=NULL,gac=NULL,De=NULL,fort=FALSE,tol=10^-10,
                                  disp=FALSE,output=FALSE,out_se=FALSE,glob=FALSE){

#        [piv,Th,Bec,gac,fv,Phi,Pp,lk,np,aic,bic] = est_multi_poly(S,yv,k,start,link,disc,difl,multi,piv,Th,bec,gac)
#
# Fit Latent Class model and some restricted versions with k classes for ordinal (NA for missing data)
# this version works with two dimensione that may contemporary affect the item responses
# 
# S    : matrix of available configurations
# yv   : frequencies of the available configurations
# k1,k2: number of latent classes for the two latent variables
# X    : matrix of covariates for the multinomial logit on the class weights
# start: type of startgine values (0 = deterministic, 1 = random, 2 = external input)
# link : type of link (0 = LC, 1 = GRM, 2 = PCM)
# disc : discriminating indices (0 = constrained, 1 = free)
# difl : difficulty levels (0 = free, 1 = additive decomposition)
# lk   : maximum log-likelihood
# piv  : weights of the latent classes
# Phi  : conditional distributions given the latent classes
# np   : number of free parameters
# bic  : Bayesian information criterion
# Th,be,ga : parameters for the model 4 (Th=Psi if start==3)
# fv   : list of items constrained
# fort : T for using fortran code for covariates, F using R code only
# tol  : relative tolerance level for convergence
# disp : TRUE for displying the likelihood evolution step by step
# output : to return additional outputs (Phi,Pp,Piv)
# out_se : to return standard errors
# glob : for global parametrization of the prior probabilities
#
# OUTPUT:
# ent  : entropy

# With k=1
	if((k1==1 & k2==1) || link==0) stop("--> for LC models use est_multi_poly")
	if(max(S,na.rm=TRUE)==1 & difl!=0){
	  warning("with binary data put difl=0\n")
	  difl = 0		
	}
# Preliminaries
# check problems with input
    cov = !is.null(X)
    if(cov){
    		X = as.matrix(X)
   		namesX = colnames(X)
   		if(glob) logit_cov = "g" else logit_cov = "m"
    }else{
		logit_cov = "m"
    }
    miss = any(is.na(S))
	ns = nrow(S); J = ncol(S)
    if(miss){
		cat("Missing data in the dataset, units and items without responses are removed\n")
		ind = which(apply(is.na(S),1,all))
		if(length(ind)>0){
			S = S[-ind,]; yv = yv[-ind]
			if(!is.null(X)) X = as.matrix(X[-ind,])
			ind = which(apply(is.na(S),2,all))
			if(length(ind)>0){
				S = S[,-ind]
				miss = any(is.na(S))
	        }
	    }
    }
    if(miss){R=1*(!is.na(S)); S[is.na(S)]=0}
	lv = apply(S,2,max)+1; lm = max(lv)
	ns = nrow(S); J = ncol(S)
	n = sum(yv)
# checks about the covariates
	if(cov){
		ncov = ncol(X)
		out = aggr_data(X,fort=fort)
		Xdis = as.matrix(out$data_dis); Xlabel = out$label; Xndis = max(out$label)
		if(glob){
			XX1dis = array(0,c(k1-1,k1-1+ncov,Xndis))
    		for(i in 1:Xndis) XX1dis[,,i] = cbind(diag(k1-1),rep(1,k1-1)%o%Xdis[i,])    		
    	}else{
	    	XX1dis = array(0,c(k1-1,(k1-1)*(ncov+1),Xndis))
    		if(k1==2) II = 1 else II = diag(k1-1)
    		for(i in 1:Xndis) XX1dis[,,i] = II%x%t(c(1,Xdis[i,]))
    	}
		if(glob){
			XX2dis = array(0,c(k2-1,k2-1+ncov,Xndis))
    		for(i in 1:Xndis) XX2dis[,,i] = cbind(diag(k2-1),rep(1,k2-1)%o%Xdis[i,])    		
    	}else{
	    	XX2dis = array(0,c(k2-1,(k2-1)*(ncov+1),Xndis))
    		if(k2==2) II = 1 else II = diag(k2-1)
    		for(i in 1:Xndis) XX2dis[,,i] = II%x%t(c(1,Xdis[i,]))
    	}
    }else{
    	ncov = 0
    	XX1dis = array(diag(k1-1),c(k1-1,k1-1,1))
    	XX2dis = array(diag(k2-1),c(k2-1,k2-1,1))
    	Xlabel = rep(1,ns)
	}
# about models
	Aggr1 = diag(k1)%x%matrix(1,1,k2)
	Aggr2 = matrix(1,1,k1)%x%diag(k2)
	if(link==1) ltype = "g" else if(link==2) ltype = "l"
	if(link == 1 || link == 2){
# items related to the first latent variable
		items1 = sort(unique(as.vector(multi1)))
		if(any(items1==0)) items1 = items1[items1>0]
		J1 = length(items1) 
		if(is.vector(multi1)) rm1 = 1 else rm1 = nrow(multi1)
# items related to the second latent variable
		items2 = sort(unique(as.vector(multi2)))
		if(any(items2==0)) items2 = items2[items2>0]
		J2 = length(items2)
		if(is.vector(multi2)) rm2 = 1 else rm2 = nrow(multi2)
# design matrix for the first latent variable
		Dem1 = matrix(0,J,rm1)
		if(rm1==1){
			Dem1 = 1
			fv1 = multi1[1]
		}else{
			for(r in 1:rm1){
				ind = multi1[r,]
				ind = ind[ind>0]
				Dem1[ind,r] = 1      
			}
			fv1 = multi1[,1]     # list of constrained items for the first latent variable
		}
		fv1e = NULL; count = 0
		for(j in 1:J){
			if(j%in%fv1) fv1e = c(fv1e,count+1)
			count = count+lv[j]-1
		}
		Dem2 = matrix(0,J,rm2)
		if(rm2==1){
			Dem2 = 1
			fv2 = multi2[1]
		}else{
			for(r in 1:rm2){
				ind = multi2[r,]
				ind = ind[ind>0]
				Dem2[ind,r] = 1      
			}
			fv2 = multi2[,1]     # list of constrained items for the second latent variable
		}
		fv2e = NULL; count = 0
		for(j in 1:J){
			if(j%in%fv2) fv2e = c(fv2e,count+1)
			count = count+lv[j]-1
		}
		fv = union(fv1,fv2)
		fve = union(fv1e,fv2e)
		rm = length(fve)
		indga1 = 1:J; indga1 = indga1[setdiff(items1,fv1)]
		indga2 = 1:J; indga2 = indga2[setdiff(items2,fv2)]
		indth1 = 1:(k1*rm1)
		indth2 = (k1*rm1+1):(k1*rm1+k2*rm2)
		if(difl==0){
			indbe = k1*rm1+k2*rm2+(1:(sum(lv-1)-rm1-rm2))
			indbec = 1:sum(lv-1); indbec = indbec[-fve]
		}else{
			indbe = k1*rm1+k2*rm2+(1:(J-rm+sum(lv[fv]-2)))
			indbec = 1:J; indbec = indbec[-fv]
		}
# abililities for each item
		abils1 = rep(0,J)
		if(rm1==1){
			abils1[multi1] = 1
		}else{
			for(h in 1:rm1){
				ind = multi1[h,]; ind = ind[ind>0]
				abils1[ind] = h
			}
		}
		abils2 = rep(0,J)
		if(rm2==1){
			abils2[multi2] = 1
		}else{
			for(h in 1:rm2){
				ind = multi2[h,]; ind = ind[ind>0]
				abils2[ind] = h
			}
		}
		abils = rep(0,J)
		for(j in 1:J) for(h in 1:rm) if(abils1[j]==abils1[fv[h]] & abils2[j]==abils2[fv[h]]) abils[j] = h
# design matrix
		if(difl==0) ZZ = array(NA,c(lm-1,k1*rm1+k2*rm2+sum(lv-1)-rm,J*k1*k2))
		if(difl==1) ZZ = array(NA,c(lm-1,k1*rm1+k2*rm2+J-rm+sum(lv[fv]-2),J*k1*k2))
		cont = 0; refitem = matrix(0,J*k1*k2,1)       # reference item of that design matrix
		for(c1 in 1:k1){
			u11 = matrix(0,1,k1); u11[c1] = 1
			for(c2 in 1:k2){
				u12 = matrix(0,1,k2); u12[c2] = 1
				for(j in 1:J){
					u21 = matrix(0,1,rm1); u21[abils1[j]] = 1
					u22 = matrix(0,1,rm2); u22[abils2[j]] = 1
					v = matrix(0,1,J); v[j] = 1
					cont = cont+1
					if(difl==0){
						Te = matrix(0,lv[j]-1,sum(lv-1))
						if(j==1) Te[,1:(lv[j]-1)] = diag(lv[j]-1)
						else Te[,sum(lv[1:(j-1)]-1)+(1:(lv[j]-1))] = diag(lv[j]-1)
						Te = matrix(Te[,-fve],lv[j]-1,dim(Te)[2]-length(fve))
					}else if(difl==1){
						Te = matrix(0,lv[j]-1,sum(lv[fv]-2))
						if(lv[j]>2){
							if(abils[j]==1) Te[,1:(lv[j]-2)] = diag(lv[j]-1)[,-1]
							else Te[,sum(lv[1:(abils[j]-1)]-2)+(1:(lv[j]-2))] = diag(lv[j]-1)[,-1]
						}
						Te = cbind(v%x%rep(1,lv[j]-1),Te)
						Te = matrix(Te[,-fv],lv[j]-1)
					}
					ZZ[1:(lv[j]-1),,cont] = cbind(rep(1,lv[j]-1)%*%(u11%x%u21),rep(1,lv[j]-1)%*%(u12%x%u22),-Te)
					refitem[cont] = j
				}
			}
		}
	}
	ZZ0 = ZZ
# inequalities for ability
	if(glob){
		II1 = diag(k1-1); II1 = cbind(0,II1)-cbind(II1,0)
		if(rm1>1) II1 = II1%x%matrix(c(1,rep(0,rm1-1)),1,rm1)
		II2 = diag(k2-1); II2 = cbind(0,II2)-cbind(II2,0)
		if(rm2>1) II2 = II2%x%matrix(c(1,rep(0,rm1-1)),1,rm1)
		II = rbind(cbind(II1,matrix(0,k1-1,dim(II2)[2])),
		           cbind(matrix(0,k2-1,dim(II1)[2]),II2))
		Dis = cbind(II,matrix(0,k1+k2-2,dim(ZZ)[2]-(k1*rm1+k2*rm2)))
	}
	else Dis = NULL
# When there is just 1 latent class
#   if k == 1,
#     piv = 1;
#     P = zeros(2,J);
#     for j in 1:J,
#       for jb = 0:1,
#         ind = which(S[,j]==jb);
#         P(jb+1,j) =  sum(yv(ind))/n;
#       end
#     end
#     Psi = P;
#     psi = ones(ns,1);
#     for j in 1:J,
#       psi = psi.*P(S[,j]+1,j);
#     end
#     lk = yv"*log(psi);
#     np = J;
#     aic = -2*lk+2*np;
#     bic = -2*lk+np*log(n);
#     bec = NULL; gac = NULL;
#     Pp = ones(ns,1);
#     Th = NULL;
#     return
#   end
# Starting values
	if(start == 0){
		if(cov){
		    de1 = de2 = NULL; Piv = matrix(1/(k1*k2),ns,k1*k2); piv = NULL				
	    }else{
			be = de1 =de2 = NULL; piv1 = rep(1,k1)/k1; piv2 = rep(1,k2)/k2
		} # latent class probabilities
		if(k1==1) grid1 = 0 else grid1 = seq(-k1,k1,2*k1/(k1-1))
		if(k2==1) grid2 = 0 else grid2 = seq(-k2,k2,2*k2/(k2-1))
		Phi = array(NA,c(lm,J,k1*k2)) # probability every response
		for(j in 1:J){
			dist = rep(0,lv[j])
			for(y in 0:(lv[j]-1)) dist[y+1] = (sum(yv[S[,j]==y])+0.5)/n
			out = matr_glob(lv[j]); Co = out$Co; Ma = out$Ma			
			eta = Co%*%log(Ma%*%dist)
			count = 0
			for(c1 in 1:k1) for(c2 in 1:k2){
				count = count+1
				Phi[1:lv[j],j,count] = inv_glob(eta+grid1[c1]+grid2[c2])$p
			}
		}
	}
	if(start == 1){
		if(cov){
			if(glob){
			    de1 = de2 = NULL; Piv = matrix(runif(ns*k1*k2),ns,k1*k2)
			    Piv = Piv*(1/rowSums(Piv))
			    piv = NULL
			}else{
				de1 = rnorm((k1-1)*(ncov+1))/rep(c(1,apply(X,2,sd)),(k1-1))
				de2 = rnorm((k2-1)*(ncov+1))/rep(c(1,apply(X,2,sd)),(k2-1))
				if(k1>1) Piv1 = prob_multi_glob(XX1dis,logit_cov,de1,Xlabel)$P 
				if(k2>1) Piv2 = prob_multi_glob(XX2dis,logit_cov,de2,Xlabel)$P
				if(k1==1) Piv = Piv2
				if(k2==1) Piv = Piv1
			    if(k1>1 & k2>1){
			    	Piv = matrix(0,ns,k1*k2)
			    	for(i in 1:ns) Piv[i,] = Piv1[i,]%x%Piv2[i,]
			    }
				piv = NULL
			}
		}else{
			piv1 = runif(k1)
			piv1 = piv1/sum(piv1)
			piv2 = runif(k2)
			piv2 = piv2/sum(piv2)
		}
		Phi = array(runif(lm*J*k1*k2),c(lm,J,k1*k2))
		for(c in 1:(k1*k2)) for(j in 1:J){
			Phi[1:lv[j],j,c] = Phi[1:lv[j],j,c]/sum(Phi[1:lv[j],j,c])
			if(lv[j]<lm) Phi[(lv[j]+1):lm,j,c] = NA
		}
		if(glob){
			if(runif(1)>0.5) for(j in 1:J){
				mPhi = (0:(lv[j]-1))%*%Phi[1:lv[j],j,]
				ind = order(mPhi)
				Phi[,j,] = Phi[,j,ind]			
			}
		}
	}
	if(start==2) de = as.vector(De)  
	# if(link==0){
		# ga = NULL
	# }else{
		if (start==0 || start==1){
			ga1 = rep(1,J1-rm1)
			ga2 = rep(1,J2-rm2)
		}else{
			ga = gac[-fv]
		}
	# }
# Compute log-likelihood
	Psi = matrix(1,ns,k1*k2) # probability observed response
	if(miss){
		for(j in 1:J) for(c in 1:(k1*k2)) Psi[,c] = Psi[,c]*(Phi[S[,j]+1,j,c]*R[,j]+(1-R[,j]))
	}else{
    	# if(fort){
            # o = .Fortran("lk_obs",J,as.integer(k),as.integer(ns),as.integer(S),as.integer(l),Phi,Psi=Psi)
            # Psi = o$Psi
        # }else{
	    	for(j in 1:J) for(c in 1:(k1*k2)) Psi[,c] = Psi[,c]*Phi[S[,j]+1,j,c]
        # }            	
	}
	if(cov){
		if(start==2){
			Piv1 = prob_multi_glob(XX1dis,ltype,de1,Xlabel)$P
			Piv2 = prob_multi_glob(XX2dis,ltype,de2,Xlabel)$P
			if(k1==1) Piv = Piv2
			if(k2==1) Piv = Piv1
		    if(k1>1 & k2>1) for(i in 1:ns) Piv[i,] = Piv1[i,]%x%Piv2[i,]
		}			
	}else{
		Piv = rep(1,ns)%o%(piv1%x%piv2)
	}
	if(k1*k2==1) Pj = Psi else Pj = Psi*Piv
	pm = rowSums(Pj)
	lk = sum(yv*log(pm))
	cat("*-------------------------------------------------------------------------------*\n")
	if(link==1 || link==2){
		cat(c("Model with multidimensional structure\n"))
		names11 = NULL
		for(j in 1:rm1){names11 = c(names11,paste("Dimension",j))}
		multi1_out = as.matrix(multi1)
		if(rm1 == 1) multi1_out = t(multi1_out)
		rownames(multi1_out) = names11
		print(multi1_out)
		names12 = NULL
		for(j in 1:rm2){names12 = c(names12,paste("Dimension",j))}
		multi2_out = as.matrix(multi2)
		if(rm2 == 1) multi2_out = t(multi2_out)
		rownames(multi2_out) = names12
		print(multi2_out)
	}
	cat(c("Link of type =                 ",link,"\n"))
	cat(c("Discrimination index =         ",disc,"\n"))
	cat(c("Constraints on the difficulty =",difl,"\n"))
	cat(c("Type of initialization =       ",start,"\n"))
	cat("*-------------------------------------------------------------------------------*\n")		
	if(disp){
		# if(link==0){
    			# cat("------------|-------------|-------------|-------------|-------------|\n");
    			# cat("  iteration |   classes   |    model    |      lk     |    lk-lko   |\n");
    			# cat("------------|-------------|-------------|-------------|-------------|\n");
    		# }else{
			if(disc==0 || (length(ga1)==0 & length(ga2)==0)){
    				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
    				cat("  iteration |   classes   |    model    |      lk     |    lk-lko   |     dis     |   min(par)  |   max(par)  |\n");
    				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
			}
			if(disc==1){
				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
				cat("  iteration |   classes   |    model    |    lk       |    lk-lko   |      dis    |   min(ga)   |   max(ga)   |   min(par)  |   max(par)  |\n");
				cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
			}
		# }
		cat(sprintf("%11g",c(0,k1*10+k2,link,lk)),"\n",sep=" | ")
	}
 	it = 0; lko = lk-10^10; dis = 0; par = NULL; dga = NULL; lkv = NULL
# Iterate until convergence
	while(((abs(lk-lko)/abs(lko)>tol) && it<10^4) || it<2){
#t0 = proc.time()
		it = it+1
		paro = par; ga1o = ga1; ga2o = ga2; piv1o = piv1; piv2o = piv2; de1o = de1; de2o = de2; lko = lk
# ---- E-step ----
		V = ((yv/pm)%o%rep(1,k1*k2))*Piv*Psi; sV = colSums(V)
#print(proc.time()-t0)
# ---- M-step ----
		YY = matrix(NA,J*k1*k2,lm)
		count = 0
		for(c in 1:(k1*k2)) for(j in 1:J){
			count = count+1
			for(y in 1:lv[j]){
				ind = (S[,j]==(y-1))
				if(miss) YY[count,y] = sum(V[ind,c]*R[ind,j]) else YY[count,y] = sum(V[ind,c])			
			}
		}
		# if(link==0){  # LC model
			# if(glob){
				# out = est_multi_glob(YY,ZZ,ltype,be=par,Dis=Dis)						
		    	# par = out$be; P = out$P		    		
				# Phi = array(t(P),c(l,J,k))
			# }else{
				# cont = 0
				# for(c in 1:k) for(j in 1:J){
					# cont = cont+1
					# Phi[,j,c] = YY[cont,]/sum(YY[cont,])
				# }
			# }
		# }else{          # other models
#print(proc.time()-t0)			
			if(disc==1){
				if(it>1 & rm<J){
					ZZ1 = ZZ2 = array(NA,c(lm-1,J,J*k1*k2))
					count = 0
					for(c1 in 1:k1) for(c2 in 1:k2) for(j in 1:J){
						count = count+1
						ZZ1[1:(lv[j]-1),,count] = 0
						ZZ1[1:(lv[j]-1),j,count] = ZZ0[1:(lv[j]-1),1:(k1*rm1),count]%*%par[1:(k1*rm1)]
						ZZ2[1:(lv[j]-1),,count] = 0
						ZZ2[1:(lv[j]-1),j,count] = ZZ0[1:(lv[j]-1),k1*rm1+1:(k2*rm2),count]%*%par[k1*rm1+1:(k2*rm2)]
					}
					ZZ1 = array(ZZ1[,indga1,],c(lm-1,length(ga1),J*k1*k2))
					ZZ2 = array(ZZ2[,indga2,],c(lm-1,length(ga2),J*k1*k2))
					if(difl==0) ind = k1*rm1+k2*rm2+1:(sum(lv-1)-rm)
					if(difl==1) ind = k1*rm1+k2*rm2+1:(J-rm+sum(lv[fv]-2))
					ZZ1Int = ZZ2Int = array(NA,c(lm-1,J*k1*k2))
					count = 0
					for(c1 in 1:k1) for(c2 in 1:k2) for(j in 1:J){
						count = count+1
						ZZ1Int[1:(lv[j]-1),count] = ga2c[j]*ZZ0[1:(lv[j]-1),k1*rm1+1:(k2*rm2),count]%*%par[k1*rm1+1:(k2*rm2)]+ZZ0[1:(lv[j]-1),ind,count]%*%par[ind]
						ZZ2Int[1:(lv[j]-1),count] = ga1c[j]*ZZ0[1:(lv[j]-1),1:(k1*rm1),count]%*%par[1:(k1*rm1)]+ZZ0[1:(lv[j]-1),ind,count]%*%par[ind]
					}
					ga1 = est_multi_glob_gen(YY,ZZ1,ltype,be=ga1,Int=ZZ1Int)$be
					ga2 = est_multi_glob_gen(YY,ZZ2,ltype,be=ga2,Int=ZZ2Int)$be
				}
				ga1c = rep(1,J); ga1c[indga1] = ga1
				ga2c = rep(1,J); ga2c[indga2] = ga2
				ZZ = ZZ0
				for(j in 1:J){
	    				ind = (refitem==j)
		    			ZZ[,1:(k1*rm1),ind] = ZZ[,1:(k1*rm1),ind]*ga1c[j]
		    			ZZ[,k1*rm1+1:(k2*rm2),ind] = ZZ[,k1*rm1+1:(k2*rm2),ind]*ga2c[j]
				}
			}
			out = est_multi_glob_gen(YY,ZZ,ltype,be=par,Dis=Dis)						
    		par = out$be; P = out$P
			Phi = array(t(P),c(lm,J,k1*k2))
		# }
#		if(any(ga1<0) | any(ga2<0)) browser()
# Update piv
		if(cov){
			if(k1>1){
				out = est_multi_glob(V%*%t(Aggr1),XX1dis,logit_cov,Xlabel,de1)
				de1 = out$be; P1dis = out$Pdis; Piv1 = out$P
			}
			if(k2>1){
				out = est_multi_glob(V%*%t(Aggr2),XX2dis,logit_cov,Xlabel,de2)
				de2 = out$be; P2dis = out$Pdis; Piv2 = out$P
			}
			if(k1==1) Piv = Piv2
			if(k2==1) Piv = Piv1
		    if(k1>1 & k2>1) for(i in 1:ns) Piv[i,] = Piv1[i,]%x%Piv2[i,]
		}else{
			piv1 = as.vector(Aggr1%*%sV)/n
			piv2 = as.vector(Aggr2%*%sV)/n
			Piv = rep(1,ns)%o%(piv1%x%piv2)
		} 
#print(proc.time()-t0)
# Compute log-likelihood
		Psi = matrix(1,ns,k1*k2)
		if(miss){
			for(j in 1:J) for(c in 1:(k1*k2)) Psi[,c] = Psi[,c]*(Phi[S[,j]+1,j,c]*R[,j]+(1-R[,j]))	
		}else{
			# if(fort){
                # o = .Fortran("lk_obs",J,as.integer(k),as.integer(ns),as.integer(S),as.integer(l),Phi,Psi=Psi)
                # Psi = o$Psi
            # }else{
    			for(j in 1:J) for(c in 1:(k1*k2)) Psi[,c] = Psi[,c]*Phi[S[,j]+1,j,c]
    		# }            	
		}
		if(k1*k2==1) Pj=Psi else Pj = Psi*Piv
        pm = rowSums(Pj)
	    lk = sum(yv*log(pm))
		if(it>1 & link>0) dis = max(c(abs(par-paro),abs(ga1-ga1o),abs(ga2-ga2o),abs(piv1-piv1o),abs(piv2-piv2o)))
		if(disp){
			if(it/10==floor(it/10)){
				# if(link==0){
					# cat(sprintf("%11g",c(it,k1*10+k2,link,lk,lk-lko)),"\n",sep=" | ")
				# }else{
					if(disc==0 || (length(ga1)==0 & length(ga2)==0)) cat(sprintf("%11g",c(it,k1*10+k2,link,lk,lk-lko,dis,min(par),max(par))),"\n",sep=" | ")
					if(disc==1) cat(sprintf("%11g",c(it,k1*10+k2,link,lk,lk-lko,dis,min(c(ga1,ga2)),max(c(ga1,ga2)),min(par),max(par))),"\n",sep=" | ")
				# }
			}
		}
		lkv = c(lkv,lk)
	}
	if(disp){
		if(it/10>floor(it/10)){
			# if(link==0){
				# cat(sprintf("%11g",c(it,k,link,lk,lk-lko)),"\n",sep=" | ")
			# }else{
				if(disc==0 || (length(ga1)==0 & length(ga2)==0)) cat(sprintf("%11g",c(it,k1*10+k2,link,lk,lk-lko,dis,min(par),max(par))),"\n",sep=" | ")
				if(disc==1) cat(sprintf("%11g",c(it,k1*10+k2,link,lk,lk-lko,dis,min(c(ga1,ga2)),max(c(ga1,ga2)),min(par),max(par))),"\n",sep=" | ")
			# }
		}
		# if(link==0){
    			# cat("------------|-------------|-------------|-------------|-------------|\n")
		# }else{
			if(disc==0 || (length(ga1)==0 & length(ga2)==0)) cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
			if(disc==1) cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
		# }
	}
# Compute number of parameters  
	np = k1*rm1+k2*rm2+disc*(J-rm)
	if(cov){
		if(glob){
			np = np+k1+k2-2+2*ncov
      	}else{
      		np = np+(k1+k2-2)*(ncov+1)
      	}
	}else{
   		np = np+k1-1+k2-1
  	}
	if(difl==0) np = np+sum(lv-1)-rm
	else if(difl==1) np = np+J-rm+sum(lv[fv]-2)
# extract parameter estimates  
	aic = -2*lk+2*np
	bic = -2*lk+np*log(n)
	th1 = par[indth1]; th2 = par[indth2]; be = par[indbe]
	if(difl==0){
		bec = rep(0,sum(lv-1)); bec[indbec] = be
		Bec = matrix(NA,J,lm-1)
		count = 0
		for(j in 1:J){
			Bec[j,(1:lv[j]-1)] = bec[count+(1:(lv[j]-1))]
			count = count+lv[j]-1			
		}
   	}else if(difl==1){
   		bec1 = rep(0,J); bec1[indbec] = be[1:(J-rm)]
   		if(rm==1){
			bec2 = rep(0,lm-1); bec2[2:(lm-1)] = be[J-rm+(1:(lm-2))]
		}else{
			bec2 = matrix(NA,lm-1,rm); bec2[1,] = 0
			count = 0
			for(h in 1:rm) if(lv[fv[h]]>2){
				bec2[2:(lv[h]-1),h] = be[J-rm+count+(1:(lv[h]-2))]
				count = count+lv[h]-2
			}
			dimnames(bec2) = list(level=1:(lm-1),dim=1:rm)
		}
		Bec = list(difficulties=bec1,cutoffs=bec2)
	}
   	ga1c = rep(1,J); ga1c[indga1] = ga1
   	ga2c = rep(1,J); ga2c[indga2] = ga2
   	Th1 = matrix(th1,rm1,k1)
   	Th2 = matrix(th2,rm2,k2)
   	names21 = NULL
   	for(c in 1:k1){names21 = c(names21,paste("Class",c))}
   	names22 = NULL
   	for(c in 1:k1){names22 = c(names22,paste("Class",c))}
   	rownames(Th1) = names11; colnames(Th1) = names21
   	rownames(Th2) = names12; colnames(Th1) = names22
	Pp = ((1./pm)%o%rep(1,k1*k2))*Piv*Psi
	Pp1 = Pp%*%t(Aggr1); Pp2 = Pp%*%t(Aggr2)
	ent = -sum(V*log(pmax(Pp,10^-100)))
	if(cov){
		if(glob){
			if(k1==1){
				De1 = NULL
			}else{
				De1 = matrix(de1,ncov+k1-1,1)
				names_cutoff = paste("cutoff",1:(k1-1),sep="")
   	   			if(is.null(namesX)){
      				namesX1 = c(names_cutoff,paste("X",1:ncov,sep=""))	
      			}else{
					namesX1 = c(names_cutoff,namesX)
	  			}
	  			rownames(De1) = namesX1
	  		}
			if(k2==1){
				De2 = NULL
			}else{
				De2 = matrix(de2,ncov+k2-1,1)
				names_cutoff = paste("cutoff",1:(k2-1),sep="")
   	   			if(is.null(namesX)){
      				namesX2 = c(names_cutoff,paste("X",1:ncov,sep=""))	
      			}else{
					namesX2 = c(names_cutoff,namesX)
	  			}
	  			rownames(De2) = namesX2
	  		}	  		
		}else{
   	   		if(is.null(namesX)){
      			namesX = c("intercept",paste("X",1:ncov,sep=""))	
      		}else{
				namesX = c("intercept",namesX)
	  		}
	  		if(k1==1) De1 = NULL
	  		else{
		 		De1 = matrix(de1,ncov+1,k1-1)
			    dimnames(De1) = list(namesX,logit=2:k1)
			}
			if(k2==1) De2 = NULL
			else{
		 		De2 = matrix(de2,ncov+1,k2-1)
			    dimnames(De2) = list(namesX,logit=2:k2)
			}
	  	}
		piv1 = colMeans(Piv1)
		piv2 = colMeans(Piv2)
	}else{
		De1 = De2 = NULL
	}
	dimnames(Phi) = list(category=0:(lm-1),item=1:J,class=1:(k1*k2))
	# if(link==0 & !glob){
		# Phi = pmax(Phi,10^-50)
  		# par = NULL
  		# for(c in 1:k) for(j in 1:J) par = c(par,log(Phi[-1,j,c]/Phi[1,j,c]))
	# }
	if(!cov){
		de1 = De1 = log(piv1[-1]/piv1[1])
		de2 = De2 = log(piv2[-1]/piv2[1])
	}
	if(out_se){
		lde1 = length(de1); lde2 = length(de2); lpar = length(par); lga = 0
		par_comp = c(de1,de2,par)
		if(disc==1){
			lga1 = length(ga1); lga2 = length(ga2)
			par_comp = c(par_comp,ga1,ga2)
		}
		if(disp){
			cat("computation of derivatives\n")
			cat(length(par_comp),"parameters\n")
		}
  		out = lk_obs_score_within(par_comp,lde1,lde2,lpar,lga1,lga2,S,R,yv,k1,k2,rm1,rm2,lv,J,fv,link,disc,indga1,indga2,
  		                         glob,refitem,miss,ltype,XX1dis,XX2dis,Xlabel,ZZ0,fort)
		scn = rep(0,length(par_comp)); Jn = NULL
		for(j in 1:length(par_comp)){
			par_comp1 = par_comp; par_comp1[j] = par_comp1[j]+10^-6
			out1 = lk_obs_score_within(par_comp1,lde1,lde2,lpar,lga1,lga2,S,R,yv,k1,k2,rm1,rm2,lv,J,fv,link,disc,indga1,indga2,
			                           glob,refitem,miss,ltype,XX1dis,XX2dis,Xlabel,ZZ0,fort)
			 scn[j] = (out1$lk-lk)*10^6	
			 Jn = cbind(Jn,(out1$sc-out$sc)*10^6)
			if(disp) if(j/10>floor(j/10)) cat(".") else cat(round(j/10))
	  		if(disp) if(j/100==floor(j/100)) cat("\n")		
  		}
  		if(disp) cat("\n")  		
  		Jn = (Jn+t(Jn))/2
  		Vn = ginv(Jn)
  		se = sqrt(abs(diag(Vn)))
  		if(k1>1) sede1 = se[1:lde1]
  		if(k2>1) sede2 = se[lde1+1:lde2]
  		separ = se[lde1+lde2+(1:lpar)]
  		if(disc==1){
  			sega1 = se[lpar+lde1+lde2+(1:lga1)]
  			sega2 = se[lpar+lde1+lde2+lga1+(1:lga2)]
  		}else{
  			sega1 = NULL; sega2 = NULL
  		}
		if(glob){
			if(k1==1) seDe1 = NULL else seDe1 = matrix(sede1,ncov+k1-1,1)
			if(k2==1) seDe2 = NULL else seDe2 = matrix(sede2,ncov+k2-1,1)
		}else{
			if(k1==1) seDe1 = NULL else seDe1 = matrix(sede1,ncov+1,k1-1)
			if(k2==1) seDe2 = NULL else seDe2 = matrix(sede2,ncov+1,k2-1)
		}	
# to check derivatives		
#		print(c(lk,out$lk,lk-out$lk))
#  		print(cbind(out$sc,scn,out$sc-scn))
	}
    out = list(piv1=piv1,piv2=piv2,Th1=Th1,Th2=Th2,Bec=Bec,ga1c=ga1c,ga2c=ga2c,fv1=fv1,fv2=fv2,
               De1=De1,De2=De2,Phi=Phi,lk=lk,np=np,aic=aic,bic=bic,ent=ent,call=match.call())
    if(output){out$Piv=Piv; out$Pp=Pp;  out$Pp=Pp1;  out$Pp=Pp2; out$lkv = lkv}
    if(out_se){out$seDe1=seDe1; out$seDe2=seDe2; out$separ=separ; out$sega1=sega1; out$sega2=sega2; out$Vn=Vn}
	class(out) = "est_multi_poly_within"
  	return(out)

}