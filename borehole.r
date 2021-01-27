#borehole: mimicing simulation code
#boreholetrue: mimicing physical experiment
#z: qualitative parameter
#zmax: a vector--the number of qualitative level for each qualitative factor


borehole<-function(x,z,zmax) {
    rw=x[,1]*(0.15-0.05)+0.05
    Tu=x[,2]*(115600-63070)+63070
    Tl=(x[,3])*(116-63.1)+63.1
    L=(x[,4])*(1680-1120)+1120
    r=(x[,5])*(50000-100)+100
    H=seq(400,410,length.out=zmax[1])[z$X1]
    Kw=seq(5000,10000,length.out=zmax[2])[z$X2]

    y=2*pi*Tu*H
    y=y/log(r/rw)
    y1=1+Tu/Tl+(2*L*Tu)/(Kw*rw^2*log(r/rw))
    y=log(y/y1)#+0.5*rnorm(nrow(x))
}


boreholetrue<-function(x) {
    rw=x[,1]*(0.15-0.05)+0.05
    Tu=x[,2]*(115600-63070)+63070
    Tl=(x[,3])*(116-63.1)+63.1
    L=(x[,4])*(1680-1120)+1120
    r=(x[,5])*(50000-100)+100
    H=401
    Kw=9000
      

    y=2*pi*Tu*H
    y=y/log(r/rw)
    y1=1+Tu/Tl+(2*L*Tu)/(Kw*rw^2*log(r/rw))
    y=log(y/y1)+rnorm(nrow(x))
}
