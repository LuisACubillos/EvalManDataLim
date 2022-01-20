
# Fase de proyeccion stochastico


nsim = 1000
nhcr = 4 #numero de reglas de control
nproy = 40 #numero de años proyectados
yrp = seq(2020 + 1, 2020 + nproy+1, 1) #años futuros
B_Bmsy = B = C = Ft = array(data = NA, dim = c(nsim, nproy+1,nhcr))
#Se seleccionan nsim valores de los parametros viables
set.seed(378)
select <- floor(runif(n = nsim,min = 1,max = 3926))
plot(OCm1$r_viable,OCm1$k_viable)

#Nota implementar transformando limites de confianza de r y k
# Xprom = (Xmin+Xmax)/2 sigma = raiz((logXprom-logXmin)«2)/2
#r.prom <-
#sigmar <- 
#k.prom <-
#sigmak <-
#corr_rk <-
  
for(l in 1:nsim){
  ll <- select[l]
  r <- OCm1$r_viable[ll]
  k <- OCm1$k_viable[ll]*1000
  Bend <- OCm1$bt_viable[ll,62]*k
  bmsy <- k/2
  fmsy <- r/2
  Yend <- OCm1$ref_ts$catch[62]
  for(j in 1:nhcr){
    for(i in 1:(nproy+1)){
      if(i == 1){
        #Proyeccion año 1
        B[l,i,j] = Bend + r*Bend*(1 - Bend/k) - Yend
        Brel = Bend/bmsy
        B_Bmsy[l,i,j]=Brel
        #APLICACION DE REGLAS
        if(j==1){
          Ft[l,i,j]=0
        }
        if(j==2){
          if(Brel <= 0.5){
            Ft[l,i,j] = 0
          }
          else{Ft[l,i,2] = fmsy}
        }
        if(j==3){
          if(Brel < 1){
            Ft[l,i,3] = Brel*fmsy
          }
          else{Ft[l,i,3] = fmsy}
        }
        if(j==4){
          if(Brel < 1){
            Ft[l,i,4] = fmsy*(Brel - 0.5)/0.5
            if(Brel <= 0.5){
              Ft[l,i,4] = 0
            }
          }
          else{Ft[l,i,4] = fmsy}
        }
        C[l,i,j] = Ft[l,i,j]*B[l,i,j]
      }else{
        # a partir del segundo año
        B[l,i,j] = max(1,B[l,i-1,j] + r*B[l,i-1,j]*(1 - B[l,i-1,j]/k) - C[l,i-1,j])
        Brel = B[l,i-1,j]/bmsy
        B_Bmsy[l,i,j]=Brel
        #APLICACION DE REGLAS
        if(j==1){
          Ft[l,i,j]=0
        }
        if(j==2){
          if(Brel <= 0.5){
            Ft[l,i,j] = 0
          }
          else{Ft[l,i,2] = fmsy}
        }
        if(j==3){
          if(Brel < 1){
            Ft[l,i,3] = Brel*fmsy
          }
          else{Ft[l,i,3] = fmsy}
        }
        if(j==4){
          if(Brel < 1){
            Ft[l,i,4] = fmsy*(Brel - 0.5)/0.5
            if(Brel <= 0.5){
              Ft[l,i,4] = 0
            }
          }
          else{Ft[l,i,4] = fmsy}
        }
        C[l,i,j] = Ft[l,i,j]*B[l,i,j]
      }
    }
    
  }
}

dim(B)
plot(yrp,B[10,,1],col=1,ty="l")
lines(yrp,B[10,,2],col=2)
lines(yrp,B[10,,3],col=3)
lines(yrp,B[10,,4],col=4)



quant <- function(x){quantile(x,c(0.025,0.5,0.975),na.rm=TRUE)}

db <- NULL
data <- B
nsim = dim(data)[1]
nyear = dim(data)[2]
nhcr = dim(data)[3]
HCR <- c("A","B","C","D")
for(iter in 1:nhcr){
  tmp_db <- cbind(mod = data[,,iter])
  tmp_db <- t(apply(t(tmp_db),1,quantile,c(0.1,0.2,0.5,0.8,0.9),na.rm=TRUE))
  db <- rbind(db,tmp_db)
}

mod <- rep(1:nhcr,each=nyear)
iyrs <- rep(c(2021:(2021+nyear-1)),nhcr)
colnames(db) <- c("Li","L1","Mediana","L2","Ls")
db_mod <- cbind(as.data.frame(db),Year=iyrs,HCR=HCR[mod])
df <- db_mod
db <- NULL

library(ggplot2)
p <- ggplot(data=df,aes(x=Year,y=Mediana,group=HCR))+
  geom_ribbon(aes(ymin = Li,ymax= Ls),fill="grey")+
  geom_ribbon(aes(ymin = L1,ymax= L2),fill="grey30")+
  geom_line()+
  facet_wrap(~HCR,ncol=4)+
  scale_y_continuous(name="Biomasa")+
  theme_bw()
p
