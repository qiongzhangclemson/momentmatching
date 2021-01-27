rm(list=ls())
require(reshape)
require(lattice)
require(nloptr)
require(glasso)
source("update.r")
source("selectx.r")
source("fullyrs.r")
source("genW.r")
source("cost.r")

set.seed(10)
K<-10#number of alternatives
rho<--0.5
IT<-100   #number of replications
N<-min(K*100,1000)     #number of samples: use a matrix to store a vector at each iteration
#for datacase 2, K has to be dividable by 3
datacase<-1 #determine how to generate data 1-multivariate normal, 2--borehole, 3--han,4--calibration borehole,5--calibration
burn<-K+3+10#warmup sample same propotion as Qu et al
#burn need to be > K+2 (b/c update formula in proposition 3)

#prior
b0<-burn-1 #number of sample used to estimate B
q0<-b0






#methods and parameters
methods<-data.frame(cbind(c(1,2,2,3,3,3,4),c(1:3,1:3,4)))
methods$method<-c("KL","Moment","MomentKL",
                      "RandomKL","RandomMoment","RandomIm","Independent")
names(methods)[1:2]<-c("selectmethod","updatemethod")
methods<-methods[c(1,2,3,7),]



cost<-NULL
bs<-matrix(0,nrow(methods),IT) # final b for each method in each replication
qs<-matrix(0,nrow(methods),IT) # final q for each method in each replication
xsall<-list() # selected x in the all replications

for(it in 1:IT) { #run all replications
	Wall<-genW(N+burn,K,datacase,rho)
    Wburn<-Wall[1:burn,]
	theta0<-apply(Wburn,2,mean)
	B0<-cov(Wburn)*(b0-K-1)
	ESigma<-array(B0/(b0-K-1),dim=c(K,K,nrow(methods)+2)) #resulted variances results
	
	
    xsall[[it]]<-matrix(0,nrow(methods),N)
    W<-Wall[-(1:burn),] #generate sample: a N*K matrix
    results<-fullyall(W,b0,q0,B0,theta0) #fully update
    mus<-results$theta #will be used for computing opportunity cost
    #mus<-apply(W,2,mean)
    ESigma[,,2]<-results$B/(b0+N-K-1)# since update all, the last b is b0+N
    
    
    data<-matrix(0,0,5) #matrix to store results from  different methods
    for(case in 1:nrow(methods)) {
        selectmethod<-methods[case,1]
        updatemethod<-methods[case,2]
        results<-fullyrs(W,b0,q0,B0,theta0,selectmethod,updatemethod)
        ESigma[,,2+case]<-results$Bs[[N]]/(results$b-K-1)
        
        bs[case,it]<-results$b
        qs[case,it]<-results$q
        xsall[[it]][case,]<-results$xs
        data<-rbind(data,results$thetas)
    }
    #for independent only, no need to devide b-K-1
    ESigma[,,nrow(methods)+2]<-results$Bs[[N]]
    
    #compute opportunity cost
    data<-merge(data,methods)
    mu<-mus[nrow(mus),]
    cost<-rbind(cost,computecost(data,mu))
    print(it)
}


cost.all<-aggregate(cost~method+step,data=cost,mean)
pdf(paste(rho,"datacase",datacase,"alternative",K,"burn",burn,"cost.pdf",sep="")) #plot cost
require(ggplot2)
myplot<-ggplot(cost.all,aes(x=step,y=cost,color=method))+geom_line(aes(linetype=method))
myplot+theme(legend.position="bottom")
dev.off()



#remove the off diagnal entry for indepnent case
ESigma[,,dim(ESigma)[3]]<-diag(diag(ESigma[,,dim(ESigma)[3]]))
cors<-function(u) {
	sds<-1/sqrt(diag(u))
	diag(sds)%*%u%*%diag(sds)
}

for(i in 1:dim(ESigma)[3]) {
    ESigma[,,i]<-cors(ESigma[,,i])
}

titles<-c("Prior","Empirical Estimation",methods$method)
mymat<-melt(ESigma)
names(mymat)[4]="Correlation"
mymat<-within(mymat,{
  X3<-factor(X3)
  levels(X3)<-titles
}
              )
pdf(paste(rho,"datacase",datacase,
          "alternative",K,"burn",burn,"cov.pdf",sep=""),
    width=9,height=6)
p<-ggplot(mymat,aes(X1,X2))+geom_tile(aes(fill=Correlation))+facet_wrap(~X3)
p+geom_text(aes(label=round(Correlation,2)),size=12/K)+xlab("")+ylab("")+scale_y_reverse()
dev.off()

selectmat<-aggregate(value~method+step,data=cost,mean)
pdf(paste(rho,"datacase",datacase,
          "alternative",K,"burn",burn,"select.pdf",sep=""))
names(selectmat)[3]="Best"
ggplot(selectmat,aes(x=step,y=method))+geom_tile(aes(fill=Best))
dev.off()
save(xsall,cost,data,bs,ESigma,file=paste(rho,"datacase",datacase,"alternative",K,"burn",burn,"data.RData",sep=""))





