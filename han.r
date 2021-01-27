han<-function(x,z) {
    #set.seed(z)
    #K<-max(z)
    #H<-matrix(runif(K*3),K,3)
    y<-cbind(-x^2,-x,1)%*%matrix(runif(3),3,1)
    y<-y+rnorm(length(x))

    y
}


hantrue<-function(x) {
    y<-cbind(-x^2,-x,1)%*%matrix(c(0.2,0.9,0.1),3,1)
    y<-y+rnorm(length(x))
    y
  
}
