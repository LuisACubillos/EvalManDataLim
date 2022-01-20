Brel <- seq(0,3,0.01)
Ft <- numeric()
Fmsy <- 0.3
# REGLA A: Ft=0
Ft = rep(0,length(Brel))
A = Ft/Fmsy
# REGLA B:Escala
for(i in 1:length(Brel)){
  if(Brel[i] <= 0.5){
    Ft[i] = 0
  }
  else{Ft[i] = Fmsy}
}
B = Ft/Fmsy
# REGLA C: Rampa 0-40
for(i in 1:length(Brel)){
  if(Brel[i] < 1){
    Ft[i] = Brel[i]*Fmsy
  }
  else{Ft[i] = Fmsy}
}
C = Ft/Fmsy
# REGLA D: Rampa 20-40
for(i in 1:length(Brel)){
  if(Brel[i] < 1){
    Ft[i] = Fmsy*(Brel[i] - 0.5)/0.5
    if(Brel[i] <= 0.5){
      Ft[i] = 0
    }
  }
  else{Ft[i] = Fmsy}
}
D = Ft/Fmsy

## GRAFICA TODAS
plot(Brel,Brel,type="n",xlab="B/Bmsy",ylab="F/Fmsy",xlim=c(0,2),ylim=c(0,2))
lines(Brel,A,lwd=2)
lines(Brel,B,col="pink",lwd=2)
lines(Brel,C,col="cyan",lwd=2)
lines(Brel,D,col="brown",lwd=2)
abline(v=0.5)
abline(v=1)
polygon(c(0,0.5,0.5,0),c(2,2,0,0),col=rgb(1,0,0,alpha=0.2), border=NA)
polygon(c(0.5,1,1,0.5),c(1,1,0,0),col=rgb(0.255,0.255,0,alpha=0.1), border=NA)
polygon(c(0.5,2.5,2.5,0.5),c(2,2,1,1),col=rgb(0.255,0.255,0,alpha=0.1), border=NA)
polygon(c(1,2.5,2.5,1),c(1,1,0,0),col=rgb(0,1,0,alpha=0.1), border=NA)
