### FASE DE PROYECCION DETERMINISTA ####

r <- OCm1$ref_pts[1,2]
k <- OCm1$ref_pts[2,2]
bmsy <- OCm1$ref_pts[5,2]
fmsy <- OCm1$ref_pts[4,2]
Bend <- OCm1$ref_ts$b[62]
Yend <- OCm1$ref_ts$catch[62]

nhcr = 4 #numero de reglas de control
nproy = 40 #numero de años proyectados
yrp = seq(2020 + 1, 2020 + nproy, 1) #años futuros
B = matrix(data = NA, nrow = nproy+1, ncol = nhcr)
C = matrix(data = NA, nrow = nproy+1, ncol = nhcr)
Ft= matrix(data = NA, nrow = nproy+1, ncol = nhcr)
for(j in 1:nhcr){
  for(i in 1:(nproy+1)){
    if(i == 1){
      #Proyeccion año 1
      B[i,j] = Bend + r*Bend*(1 - Bend/k) - Yend
      Brel = Bend/bmsy
      #APLICACION DE REGLAS
      if(j==1){
        Ft[i,j]=0
      }
      if(j==2){
        if(Brel <= 0.5){
          Ft[i,j] = 0
        }
        else{Ft[i,2] = fmsy}
      }
      if(j==3){
        if(Brel < 1){
          Ft[i,3] = Brel*fmsy
        }
        else{Ft[i,3] = fmsy}
      }
      if(j==4){
        if(Brel < 1){
          Ft[i,4] = fmsy*(Brel - 0.5)/0.5
          if(Brel <= 0.5){
            Ft[i,4] = 0
          }
        }
        else{Ft[i,4] = Fmsy}
      }
      C[i,j] = Ft[i,j]*B[i,j]
    }else{
      # a partir del segundo año
      B[i,j] = max(1,B[i-1,j] + r*B[i-1,j]*(1 - B[i-1,j]/k) - C[i-1,j])
      Brel = B[i-1,j]/bmsy
      #APLICACION DE REGLAS
      if(j==1){
        Ft[i,j]=0
      }
      if(j==2){
        if(Brel <= 0.5){
          Ft[i,j] = 0
        }
        else{Ft[i,2] = fmsy}
      }
      if(j==3){
        if(Brel < 1){
          Ft[i,3] = Brel*fmsy
        }
        else{Ft[i,3] = fmsy}
      }
      if(j==4){
        if(Brel < 1){
          Ft[i,4] = fmsy*(Brel - 0.5)/0.5
          if(Brel <= 0.5){
            Ft[i,4] = 0
          }
        }
        else{Ft[i,4] = Fmsy}
      }
      C[i,j] = Ft[i,j]*B[i,j]
    }
  }
}

B1 <- c(OCm1$ref_ts$b,B[1:nproy,1])
B2 <- c(OCm1$ref_ts$b,B[1:nproy,2])
B3 <- c(OCm1$ref_ts$b,B[1:nproy,3])
B4 <- c(OCm1$ref_ts$b,B[1:nproy,4])
Year <- c(OCm1$ref_ts$year,yrp)
Biom <- data.frame(Year,B1,B2,B3,B4)

plot(Year,Biom$B1,ty="l",ylim=c(0,50000),ylab="Biomasa")
lines(Year,Biom$B2,col=3)
lines(Year,Biom$B3,col=4)
lines(Year,Biom$B4,col=5)
abline(h=bmsy,lty=3)
abline(h=0.5*bmsy,lty=3,col="red")
legend("topleft",c("F=0","Escala","F0-40","F20-40"),col=c(1,3,4,5),lty=c(1,1,1,1),cex=0.6)


plot(Year,Biom$B1/bmsy,ty="b",cex=0.3,ylim=c(0,2),ylab="B/Bmsy")
lines(Year,Biom$B2/bmsy,col=3)
lines(Year,Biom$B3/bmsy,col=4)
lines(Year,Biom$B4/bmsy,col=5)
abline(h=1,lty=3)
abline(h=0.5,lty=3,col="red")
legend("topleft",c("F=0","Escala","F0-40","F20-40"),col=c(1,3,4,5),lty=c(1,1,1,1),cex=0.6)