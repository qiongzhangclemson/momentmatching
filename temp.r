temp<-function(x,k) {
 #function to generate temperature of a solid sphere  samples
 #n: number of unique data points
 #ms: number of replications of each point
 #require("lhs")
 #u<-randomLHS(n,9)
 #x<-u
 n<-nrow(x)
 x[,2]<-x[,2]*2000
 x[,3]<-x[,3]*20+250
 x[,4]<-x[,4]*70-100
 x[,5]<-x[,5]*30+180
 x[,6]<-x[,6]*70+30
 x[,7]<-x[,7]*0.1+0.05
 x[,8]<-x[,8]*200+300
 x[,9]<-x[,9]*2000+7000
 pi1<-x[,5]*x[,7]/x[,6]
 pi2<-x[,6]*x[,2]/(x[,8]*x[,9]*x[,7]^2)
 y0<-rep(0,n)
 for(i in 1:k) {
   eta<-c()
   for(pit in pi1) {
      etafun<-function(u) 1-u/tan(u)-pit
      eta<-c(eta,uniroot(etafun, lower =(i-1)*pi+0.000001,
                         upper = i*pi-0.00001)$root)
    }
    yt<-exp(-eta^2*pi2)*sin(eta*x[,1])/(eta*x[,1])
    yt<-yt*4*(sin(eta)-eta*cos(eta))/(2*eta-sin(2*eta))
    y0<-y0+yt
 }
 y<-x[,3]+x[,4]*y0
 y+rnorm(nrow(x))
 #list(y=y,x=u[,1:5])

}




