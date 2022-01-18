##############################################################
#  LAB 1 GUIA PRACTICA - ONLY CATCH MODELS
#  Curso Evaluación y Manejo de Pesquerías Limitadas en Datos
#  Prof. Luis Cubillos - Universidad de Concepción
#  Enero 2022
#  Modelos de Producción Dinámicos
#  B_{t+1} = B_{t} + r B_{t} * (1-B_{t}/k) - C_{t}
##############################################################


# Package a utilizar ------------------------------------------------------
#devtools::install_github("cfree14/datalimited2")
library(datalimited2)
library(ggplot2)
#parametros gráficos
oldpar <- par()

# Datos -------------------------------------------------------------------

reineta <- read.csv("Datos/Capturas_Brama_australis.csv")
head(reineta)
colnames(reineta) <- c("yr","ct")
  

# Grafico de capturas  ----------------------------------------------------

yt <- reineta$ct
captura_reineta <- ts(yt,start = 1994,frequency = 1)
plot(captura_reineta,ty="b",pch=19,
     ylim=c(0,max(reineta$ct)*1.2),las=1,
     ylab="Captura (toneledas)",xlab="Años")
# vaersión ggplot
fig1 <- ggplot(data=reineta, aes(x = yr, y = ct))+
  geom_point()+
  geom_line()+
  ylab("Capturas (toneladas)")+
  xlab("Años")+
  theme_bw()
fig1


# Modelo 1: CMSY2 - Froese et al. (2017) ----------------------------------
catch = reineta$ct
yr = reineta$yr

m1 = cmsy2(year=yr, catch=catch, resilience="High")

# Resumen
summary(m1)

#Resumen parámetros
m1$ref_pts
r <- m1$ref_pts$est[1]
k <- m1$ref_pts$est[2]
msy <- m1$ref_pts$est[3]
fmsy <- m1$ref_pts$est[4]
bmsy <- m1$ref_pts$est[5]
r*k/4
r/2
k/2

#Resumen series de tiempo
head(m1$ref_ts)

# Diagnóstico del ajuste
plot_dlm(m1)
par <- oldpar

# Gráfico de Biomasa con límites de confianza
m1.Bt = m1$ref_ts$b
m1.Bt_l = m1$ref_ts$b_lo   
m1.Bt_h = m1$ref_ts$b_hi
plot(yr, m1.Bt, type="b", bty="n", pch=19,xlab="",ylab="Biomasa",ylim=c(0,max(m1.Bt_h)))
lines(yr, m1.Bt_l, lty=2)
lines(yr, m1.Bt_h, lty=2)

# versión ggplot
df <- data.frame(yr=yr,Bt=m1.Bt,Li=m1.Bt_l,Ls=m1.Bt_h)
fig2 <- ggplot(data=df,aes(x=yr,y=Bt))+
  geom_ribbon(aes(ymin=Li,ymax=Ls),fill="grey70")+
  geom_line()+
  geom_hline(yintercept = bmsy,linetype="dashed")+
  geom_hline(yintercept = bmsy/2,linetype="dashed")+
  scale_y_continuous(name = "Biomasa",limits = c(0,150000))+
  theme_bw()
fig2


# Modelo 2: OCOM Modelo de producción de Zhou et al. (2017a) ---------------

m2 = ocom(year = yr, catch = catch, m = 0.35)

#Resumen
summary(m2)
plot_dlm(m2)
par <- oldpar

#Resumen parámetros
m2$ref_pts
bmsy <- m2$ref_pts$q0.5[4]
#Resumen series de tiempo
m2$ref_ts

# Gráfico de Biomasa con límites de confianza
m2.Bt = m2$ref_ts$b
m2.Bt_l = m2$ref_ts$b_lo   
m2.Bt_h = m2$ref_ts$b_hi
plot(yr, m2.Bt, type="b", bty="n", pch=19,xlab="",ylab="Biomasa",ylim=c(0,max(m2.Bt_h)))
lines(yr, m2.Bt_l, lty=2)
lines(yr, m2.Bt_h, lty=2)

df <- data.frame(yr=yr,Bt=m2.Bt,Li=m2.Bt_l,Ls=m2.Bt_h)
fig3 <- ggplot(data=df,aes(x=yr,y=Bt))+
  geom_ribbon(aes(ymin=Li,ymax=Ls),fill="grey70")+
  geom_line()+
  geom_hline(yintercept = bmsy,linetype="dashed")+
  geom_hline(yintercept = bmsy/2,linetype="dashed")+
  scale_y_continuous(name = "Biomasa",limits = c(0,300000))+
  theme_bw()
fig3


# Modelo 3 : Zbrt de Zhou et al. (2017b) ----------------------------------

m3 <- zbrt(year=yr, catch=catch)

#Resumen
summary(m3)
plot_dlm(m3)

#Resumen series de tiempo
m3$ts


