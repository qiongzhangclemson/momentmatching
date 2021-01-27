fullyrs<-function(W,b0,q0,B0,theta0,selectmethod,updatemethod) {
      #update one alternative a time
      b=b0
      q=q0
      theta=theta0
      B=B0
      #thetas:matrix to collect updated mean in each step
      thetas=matrix(0,0,K)
      #xs:vector to collect selected alternative in each step
      xs<-c()
      #dbs: vector to collect increasing step of b in each step
      dbs<-c()
      #Bs: list to collect B in each step
      Bs<-list()
      for(i in 1:nrow(W)) {
          x<-selectx(q,b,theta,B,B0/(b0-K-1),selectmethod)
          xs<-c(xs,x)
          Y=W[i,]
          results=update(x,q,b,theta,B,B0/(b0-K-1),Y,updatemethod)
          q=results$q
          b=results$b
          theta=results$theta
          dbs<-c(dbs,results$db)
          B<-results$B
          Bs[[i]]<-B
          thetas <- rbind(thetas,theta)
      }
      rownames(thetas)=1:nrow(W)
      thetas=melt(thetas)
      names(thetas)<-c("step","alternative","value")
      thetas$selectmethod<-selectmethod
      thetas$updatemethod<-updatemethod
      return(list(thetas=thetas,Bs=Bs,dbs=dbs,xs=xs,b=b,q=q))
}


fullyall<-function(W,b0,q0,B0,theta0) {#update all
    b=b0
    q=q0
    B=B0
    K<-length(theta0)
    theta<-matrix(theta0,1,K)
    for(i in 1:nrow(W)) {
       thetanew<-(q*theta[i,]+W[i,])/(q+1)
       B<-B+q/(q+1)*outer(theta[i,]-W[i,],theta[i,]-W[i,],"*")
       b<-b+1
       q<-q+1
       theta<-rbind(theta,thetanew)
    }
    return(list(theta=theta,B=B))
}
