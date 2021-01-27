#x:index of the selected model
#Y:entire observation
#q,b,theta,B given parameters
#method: 1-4
#1: Qu update
#2: moment matching
#3: Improved KL
#4: Independent

update<-function(x,q,b,theta,B,B0,Y,method) {
     if(method==1) {#Qu update
         y=Y[x]
         K=length(theta)
         db=1/K
         qnew=q+db
         bnew=b+db
         #(18) in Qu
         thetanew<-(y-theta[x])/((q*bnew/(bnew-K+1)+1)*B[x,x])
         thetanew<-thetanew*B[,x]+theta
         temp=(q*(y-theta[x])^2)/(q*bnew/(bnew-K+1)+1)-B[x,x]/b
         Bnew=bnew/b*B+bnew/(b+1)*temp*outer(B[,x],B[,x],"*")/(B[x,x]^2)
     }
     if(method==2) {#moment matching
         y=Y[x]
         K=length(theta)
         db<-1/K
         bnew=b+db
         qnew=q+db
         thetanew=theta+(y-theta[x])/(B[x,x]*(q+1))*B[,x]
         tq<-(1+q*(y-theta[x])^2/((q+1)*B[x,x]))
         Bnew<-B
         Bnew[-x,x]=qnew*(bnew-K-1)*tq*B[-x,x]/((q+1)*(b-K))
         Bnew[x,-x]=Bnew[-x,x]
         Bnew[x,x]=qnew*(bnew-K-1)*tq*B[x,x]/((q+1)*(b-K))
         Bc<-B[-x,-x]-outer(B[-x,x],B[-x,x],"*")/B[x,x]
         Bnew[-x,-x]=tq*(Bc/(b-K)+outer(B[-x,x],B[-x,x],"*")/B[x,x])/(q+1)
         Bnew[-x,-x]=Bnew[-x,-x]+Bc/q     
         Bnew[-x,-x]=Bnew[-x,-x]*qnew*(bnew-K-1)/(b-K)
      }
     if(method==3) {#new update
         db<-1/K
         y=Y[x]
         K=length(theta)
         #bnew, qnew thetanew same as method 2
         bnew=b+db
         qnew=q+db
         thetanew=theta+(y-theta[x])/(B[x,x]*(q+1))*B[,x]
         Bnew=B
         Bc<-B[-x,-x]-outer(B[-x,x],B[-x,x],"*")/B[x,x]
         Bnew[x,x]<-qnew*(bnew-K+1)*(B[x,x]+q*(y-theta[x])^2/(q+1))/((b+1)*(q+1))
         Bnew[-x,x]=Bnew[x,x]*B[-x,x]/B[x,x]
         Bnew[x,-x]=Bnew[x,x]*B[x,-x]/B[x,x]
         Bnew[-x,-x]=bnew*qnew*Bc/(b*q)
         Bnew[-x,-x]=Bnew[-x,-x]+outer(Bnew[-x,x],Bnew[-x,x],"*")/Bnew[x,x]
         #Bnew<-glasso(Bnew,rho=0.001)$w
       }
       if(method==4) {#independent update
         B0<-B0/(b0-K-1)
         #arbitrarily give two numbers
         bnew<-1
         qnew<-1
         Bnew<-B
         thetanew<-theta
         y<-Y[x]
         Bnew[x,x]<-1/(1/B[x,x]+1/B0[x,x])
         thetanew[x]<-theta[x]/B[x,x]+y/B0[x,x]
         thetanew[x]<-thetanew[x]*Bnew[x,x]
        }  

     list(b=bnew,q=qnew,theta=thetanew,B=Bnew)

 }
