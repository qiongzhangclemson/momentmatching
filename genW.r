#input: number of sample, number of alternative, example case (1,2,3)
#output: an N*K matrix, containing the samples

genW<-function(N,K,case,rho) {
    if(case==1) {#W generated from multivariate normal
      thetat<-seq(0,1,len=K)             #mean of multivariate normal
      Sigma<-toeplitz(rho^seq(0,K-1,1))         #covariance matrix
      W<-matrix(rnorm(N*K),N,K)%*%chol(Sigma)
      W<-W+matrix(thetat,N,K,byrow=T)      
    }
    if(case==2) {#W generated from borehole
      #The first qualitative factor is not effective
      source("slhd.r")
      source("borehole.r")
      require("AlgDesign")
      W<-matrix(0,N,0)
      p<-5 #number of quantitative variable,i.e. dim(x)
      m<-3 #number of data points for each simulation sample
      zmax<-c(3,K/3)
      zs<-gen.factorial(zmax,length(zmax),center=F)
      for(k in 1:K) {
            x<-slhd(m,N,p)$D
            #the qualitative levels of the k-th alternative           
            z<-zs[k,]
            ws<-apply(x,3,function(xmat) borehole(xmat,z,zmax))
            W<-cbind(W,apply(ws,2,mean))
      }
    }
    if(case==3) {#W generated from simple function
      source("han.r")
      source("slhd.r")
      W<-matrix(0,N,0)
      p<-1
      m<-3
      for(k in 1:K) {
      	 x<-slhd(m,N,p)$D
         ws<-apply(x,3,han,z=k)
         W<-cbind(W,apply(ws,2,mean))
      }  
    }  
    if(case==4) {#W generated from borehole
      #The first qualitative factor is not effective
      source("slhd.r")
      source("borehole.r")
      require("AlgDesign")
      W<-matrix(0,N,0)
      p<-5 #number of quantitative variable,i.e. dim(x)
      m<-20 #number of data points for each simulation sample
      zmax<-c(3,K/3)
      zs<-gen.factorial(zmax,length(zmax),center=F)
      for(k in 1:K) {
            #the qualitative levels of the k-th alternative   
            x<-slhd(m,N,p)$D
      		wtrue<-apply(x,3,boreholetrue)        
            z<-zs[k,]
            ws<-apply(x,3,function(xmat) borehole(xmat,z,zmax))
			W<-cbind(W,-apply((wtrue-ws)^2,2,mean))
      }
    }
    if(case==5) {#W generated from simple function
      source("han.r")
      source("slhd.r")
      W<-matrix(0,N,0)
      p<-1
      m<-3
      for(k in 1:K) {
      	 x<-slhd(m,N,p)$D
         ws<-apply(x,3,han,z=k)
         wtrue<-apply(x,3,hantrue)
          W<-cbind(W,-apply((wtrue-ws)^2,2,mean))
      }  
    }  
    return(W)
}
