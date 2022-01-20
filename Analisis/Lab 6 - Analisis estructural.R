##############################################################
#  LAB 6 GUIA PRACTICA - Análisis estructural de capturas
#  Curso Evaluación y Manejo de Pesquerías Limitadas en Datos
#  Prof. Luis Cubillos - Universidad de Concepción
#  Enero 2022
##############################################################

# packages ----------------------------------------------------------------

library(strucchange)
library(stepR)
library(tseries)

# Datos de captura --------------------------------------------------------

reineta <- data.frame(
  yr = 1994:2020,
  ct = c(1186, 3930, 5585, 5998, 6332, 6828, 8159, 15156, 4429, 2645, 3764, 12707, 2517, 3743, 6160, 15199, 16977, 28814, 23079, 11955, 35975, 34218, 27586, 25267, 28175, 44288, 38109))
yt <- reineta$ct
captura_reineta <- ts(yt,start = 1994,frequency = 1)
plot(captura_reineta,ty="b",pch=19,ylim=c(0,max(reineta$ct)*1.2),las=1,ylab="Captura (toneledas)",xlab="Años")


# Clasificación del estatus -----------------------------------------------

ct <- reineta$ct
yrel <- ct/max(ct)
yrs <- 1994:2020
t_max = (yrs)[which(yrel==max(yrel))]

plot(yrs,yrel,pch=19,xlab="Años",ylab="Yt/Ymax",las=1)
abline(h=0.5,lty=2,col="grey30")
abline(h=0.1,lty=3,col="grey30")
# Si en los año anteriores a yt/ymax, entonces
# Subdesarrolo si yt<0.1
polygon(x=c(1994,1994,1995.5,1995.5),y=c(0,1,1,0),col="grey")
polygon(x=c(1995.5,1995.5,2010.5,2010.5),y=c(0,1,1,0),col="green")
polygon(x=c(2010.5,2010.5,2020,2020),y=c(0,1,1,0),col="yellow")
polygon(x=c(2010.5,2010.5,2020,2020),y=c(0.75,1,1,0.75),col="yellow")
text(2015,0.9,"¿MSY?")
text(2015,0.3,"Explotación plena")
text(2005,0.9,"Subexplotación")
points(yrs,yrel,pch=19)
abline(h=0.5,lty=2,col="grey30")
abline(h=0.1,lty=3,col="grey30")



# Evalua estacionaridad ---------------------------------------------------

yt <- log(reineta$ct)
tsy <- ts(yt,start=c(1994,1),frequency = 1)
plot(tsy,ylab="log-captura")

# Test de Dickey-Fuller
adf.test(tsy,alternative="stationary")


# Análisis estructural ----------------------------------------------------

#Test F de Chow
m0 <- Fstats(tsy~1)
sctest(m0)



# Puntos de quiebre -------------------------------------------------------

bp <- breakpoints(tsy~1)
plot(bp)
summary(bp)


# bloques de años ---------------------------------------------------------

plot(tsy)
lines(bp)
iconf <- confint(bp)
lines(iconf)

# Test de segmentación
m0b <- stepFit(y=yt,x=yrs,alpha=0.01,jumpint = TRUE,confband = TRUE)
plot(yrs,yt,pch=16,col="grey30",las=1,ylab="log(Rt)",xlab="Años")
lines(m0b,lwd=3,col="darkgrey")


# Grafico final (un salto) ------------------------------------------------

nyr = length(yrs) 
salto1 <- 15
r1 = mean(yt[1:salto1])
sd1 = sd(yt[1:salto1])
r2 = mean(yt[(salto1+1):nyr])
sd2 = sd(yt[(salto1+1):nyr])
xm <- yrs
ym <- c(rep(r1,salto1),rep(r2,(nyr-salto1)))
lsup1 <- ym[1:salto1]+sd1
linf1 <- ym[1:salto1]-sd1
lsup2 <- ym[(salto1+1):nyr]+sd2
linf2 <- ym[(salto1+1):nyr]-sd2
plot(xm,yt,ty="n",las=1,ylab="log-Captura",xlab="Años",pch=19,cex=0.8)
polygon(c(xm[1:salto1],rev(xm[1:salto1])),c(linf1,rev(lsup1)),col="gray80",border=NA)
polygon(c(xm[(salto1+1):nyr],rev(xm[(salto1+1):nyr])),c(linf2,rev(lsup2)),col="gray80",border=NA)
points(xm,yt,pch=19,cex=0.8)
lines(xm[1:salto1],ym[1:salto1],lwd=1.5)
lines(xm[(salto1+1):nyr],ym[(salto1+1):nyr],lwd=1.5)


###### a) Rendimiento máximo sostenible???
RMS <- exp(r2)

# Calcule F40
# Calcule YPR(F40)
# Calcule Reclutamiento en MSY : R = RMS/YPRf40
# Calcule la biomasa desovante inexplotada: B0 = SPRo*R

##### b) Utilice r2 como una regla de control de captura
# Pasar a Lab 7







