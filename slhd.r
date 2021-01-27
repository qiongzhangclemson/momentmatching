 ############################################################################
 #  sliceLHD: Function to generate sliced Latin Hypercube Design
 #  Arguments:
 #        m: the sample size of the small LHD
 #        n: the sample size of the large LHD, where n = mt;
 #         n is not to be provided
 #        t: the number of slices
 #        q: the number of factors
 #    Value: a n-by-q LHD with t slices forming a small LHD of size m
 #
 #    Last Update: 06/07/2010, Youngdeok Hwang
 ############################################################################

  slhd<-function(m, t, q){
      n<-m*t
      P<-matrix(0, m, t)
      Q<-matrix(0, m, t)
      xmat<-array(0,c(m,q,t))
          Pout<-array(0,c(m,q,t))
      for (l in 1:q){
          for (j in 1:t) {
             P[,j]<-sample(1:m)
          }
          for (i in 1:m) {
             Q[i,]<-sample(1:t)
          }
          newP<-P
          for (j in 1:t){
             for (i in 1:m){
                for (k in 1:m){
                   if (P[i,j]==k) { newP[i,j]<-(k-1)*t + Q[k,j] }
                }
             }
          }
          U<-runif(n, 0, 1/n)
          xmat[,l,]<-(newP-1)/n+U
                  Pout[,l,]<-newP
       }
           #xmat1<-xmat[,,1]
           #for(i in 2:t) {
           #   xmat1<-rbind(xmat1,xmat[,,i])
           #}
           #xmat<-xmat1
           out <- NULL
           out$D <- xmat
           out$P <- Pout
           out
  }
