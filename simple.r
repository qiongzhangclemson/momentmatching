simple<-function(x,z) {

   #y=cos(6.8*pi*x/2)*(z==1)-cos(7*pi*x/2)*(z==2)+cos(7.2*pi*x/2)*(z==3)+
   y=cos((6+z*0.2)*pi*x[,1]/2)*(-1)^z

   y

}