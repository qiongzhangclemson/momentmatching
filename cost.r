computecost<-function(data,mu) {
   #mu: the final empirical results
   #data: results in every step for every method
   cc<-max(mu)
   cost<-aggregate(value~method+step,data=data,which.max)
   cost$cost<-cc-mu[cost$value]
   return(cost)
}
