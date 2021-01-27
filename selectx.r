#method valued at 1--4
#1:Qu's select
#2:ZS's select
#3:random selection
#4:selection for independent case
selectx<-function(q,b,theta,B,B0,method) {
    K<-length(theta)
    if(method==4)  {#selection for independent case
       #need to recover the real variance value
       sigma2new<-1/(1/diag(B)+1/diag(B0))
       #the first term in vn(k)
       tsigma<-sqrt(diag(B)-sigma2new)
       #create a matrix with i-th row be the theta vector without the i-th entry 
       jacktheta<-matrix(theta,K,K+1,byrow=TRUE)[,2:K]
       psi<--abs(theta-apply(jacktheta,1,max))/tsigma
       fpsi<-dnorm(psi)+psi*pnorm(psi)
       #value of information
       vkg<-tsigma*fpsi
       xopt<-which.max(vkg)
    }
    if(method==3) { #random select
        xopt<-sample(K,1)
    }
    if(method==1|method==2) {
	df=b-K+1
	qk=(q+1)/(q*(b-K+1))
        if(method==1) { #qu's select
            #b^(n+1)=b^n+1/K
            S0=sqrt(qk)/(q*(b+1/K)/(b+1/K-K+1)+1)
        }
        if(method==2) { #proposed select
            S0=1/sqrt((q+1)*q*(b-K+1))
        }
        alpha=theta
        #vs:a vector collect value of informations
        vs=c()
        for(x in 1:K) {
            S<-S0/sqrt(B[x,x])*B[,x]
            beta<-S
            stopifnot(sum(S==S[1])<length(S))             
            re<-subway(alpha,beta)
            cc<-re$cc
            idx<-re$idx
            bb<-re$b[idx]
            bdeltavalue<-bb[-1]-bb[-c(length(bb))]
            expectedvalues<-dt(abs(cc),df)*(df+cc^2)/(df-1)
            expectedvalues<-expectedvalues-abs(cc)*(1-pt(abs(cc),df))
            vs=c(vs,sum(bdeltavalue*expectedvalues))
        }
		xopt=which.max(vs)
    }
    xopt
}




#input a & b: powell and ilya book p 101 after 5.19
#output:
#cc: ci values in (24) Qu et al
#b: bi values in (24) Qu et al
#idx: nondominated bi indexes
subway<-function(a,b) {
    K<-length(a)
    ab<-data.frame(cbind(a,b))
    ab<-ab[!duplicated(b),]
    K<-nrow(ab)
    ab<-ab[order(ab$b),]
    cc<-(ab[1,1]-ab[2,1])/(ab[2,2]-ab[1,2])
    idx<-c(1,2)
    if(K>2) {
     for(i in 2:(K-1)) {
       flag=0
       while(flag==0) {
           L=length(idx)
           M=L-1
           cct<-(ab[idx[L],1]-ab[i+1,1])/(ab[i+1,2]-ab[idx[L],2])
           if(cct>cc[M]) {
               flag=1
               cc=c(cc,cct)
               idx=c(idx,i+1)
           } else {
               idx<-idx[-L]
               cc<-cc[-M]
               if(length(cc)==0) {
                   flag=1
                   cc<-cct
                   idx<-c(idx,i+1)
               }
           }

       }
     }
    }
    return(list(cc=cc,idx=idx,b=ab$b))
}
