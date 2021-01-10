rm(list = ls())
#把库引入
#LIBRARIES
library(ggplot2)
library(gridExtra)    #grid.arrange()
library(reshape2)     #melt()
library(directlabels)
library(msm)          #ptnorm()
library(signal)       #interp1()
library(RColorBrewer)
library(splancs)  #inpip()


#SETUP
wd <- getwd()
outd <- "outputs"
if(!file.exists(outd)) dir.create(outd)
figd <- "figures"
if(!file.exists(figd)) dir.create(figd)
##以上都是准备工作##


##1. 经济部分

rho.T <- function(T) return(1030-0.1625*T-0.00269*T^2) #kg/m3
h.T <- function(T, cw) return((T+273)*cw)   #J/kg, T in degree Celsius
mu.T <- function(T) return(343.18e-7*10^(247.8/(T.inj+140))) #viscosity McDermott et al. (2006)
##计算V.opt
calc_dV.opt <- function(z, I.well, I.res, N.res, eta.el, rho.inj, rho.prod, cw, T0, T.gradient, T.inj, g) {
  a <- -I.well
  b <- -I.res/N.res
  c <- eta.el*rho.inj*cw*(T0+z*1e-3*T.gradient-T.inj) + z*g*(rho.inj-rho.prod)
  
  a1 <- 2.75*a   #  a1 <- 1/4*2.75*a
  b1 <- 2*b
  c1 <- c
  V0.p <- (-b1+sqrt(b1^2-4*a1*c1))/(2*a1)
  V0.m <- (-b1-sqrt(b1^2-4*a1*c1))/(2*a1)
  
  if(V0.p > 0) V0 <- V0.p else V0 <- V0.m
  a2 <- 1/2 * 2.75*1.75*0.75*a*V0^(-0.25)
  b2 <- 2.75*1.75*a*V0^(0.75) + 2*b
  c2 <- 2.75*a*V0^(1.75) + 2*b*V0 + c
  DV.p <- (-b2+sqrt(b2^2-4*a2*c2))/(2*a2)
  DV.m <- (-b2-sqrt(b2^2-4*a2*c2))/(2*a2)
  if(DV.p > 0) DV <- DV.p else DV <- DV.m
  
  Vopt <- V0+DV
  return(Vopt)
}
##计算V.cap
calc_dV.max <- function(dV.opt, sigma, z, I.res, g, Temp) {
  nu <- 0.25
  Sv <- 2500*g*z
  P.hydro <- rho.T(Temp)*g*z
  ratio <- nu/(1-nu)
  #Smin-P.hydro = nu/(1-nu) (Sv-P.hydro)
  Smin <- (1-ratio)*P.hydro+ratio*Sv
  
  dV.max <- (Smin*(1-3*sigma)-P.hydro)/I.res
  if(dV.max > dV.opt) dV.max <- dV.opt
  return(dV.max)
}

W0 <- function(y) {
  x <- y - y^2 + 3/2*y^3 - 8/3*y^4 + 125/24*y^5
  return(x)
}


##输入参数
T0 <- 15                      #deg C
cw <- 4180           #J/(kg.K)
g <- 9.81            #m/s2
I.res <- 0.2*1e9

I.w0 <- 1869*(1e3)^1.75     #Pa/(m3/s)^1.75   to check (Basel)
a.inj <- 2
Z0 <- 2500           #m
Dw <- 0.25           #m
sigma <- 0   #0.05   #0
U <- 8                     #W/m/K
Ts <- 50                 #supply temp
DTs <- Ts-T0
n.hr <- 8000    #hrs/yr
C.frac <- 1e6												#USD/well
C1 <- 151      #EUR/m
C2 <- 1378     #EUR/m2
D.DH <- 0.5   #m

type <- "triplet"

z <- seq(4, 9, 1)*1e3         #depth m
L <- seq(0, 80, 1)*1e3        #m

T.gradient.min <- 30       #deg C
T.gradient.max <- 40       #deg C
T.inj.min <- 60                   #deg C
T.inj.max <- 75                   #deg C
Dt.EGS.min <- 20    #yr
Dt.EGS.max <- 30    #yr
Dt.well.min <- 5   #yr
Dt.well.max <- 30   #yr


nz <- length(z) #6
nL <- length(L) #81
nrdm <- 1e4  #10000

##6*10000
E.el <- array(NA, dim=c(nz,nrdm))
Dt.EGS <- array(NA, dim=c(nz,nrdm))
C.well <- array(NA, dim=c(nz,nrdm))
C.EGS <- array(NA, dim=c(nz,nrdm))
E.EGS <- array(NA, dim=c(nz,nrdm))
P.EGS <- array(NA, dim=c(nz,nrdm))
#6*10000*81
E.DH <- array(NA, dim=c(nz,nrdm,nL))
C.DH <- array(NA, dim=c(nz,nrdm,nL))
En.DH <- array(NA, dim=c(nz,nrdm,nL))
E.comb <- array(NA, dim=c(nz,nrdm,nL))
C.comb <- array(NA, dim=c(nz,nrdm,nL))
P.comb <- array(NA, dim=c(nz,nrdm,nL))

for(i in 1:nz){ #6
    for(rdm in 1:nrdm) { #10000
  
        if(type == "doublet") {N.res <- 1; N.inj <- 1; N.prod <- 1; N.well <- 2}
        if(type == "triplet") {N.res <- 2; N.inj <- 2; N.prod <- 1; N.well <- 3}
        if(type == "3stageDoublet") {N.res <- 3; N.inj <- 1; N.prod <- 1; N.well <- 2}
        T.gradient <- runif(1, min=T.gradient.min, max=T.gradient.max)
        T.inj <- runif(1, min=T.inj.min, max=T.inj.max)
        
        T.prod <- T0+T.gradient*z[i]*1e-3
        h.inj <- h.T(T.inj, cw)
        h.prod <- h.T(T.prod, cw)
        rho.inj <- rho.T(T.inj)
        rho.prod <- rho.T(T.prod)
        mu.inj <- mu.T(T.inj)
        mu.prod <- mu.T(T.prod)
        
        dI.inj <- 0.241 / Dw^4.75 * mu.inj^0.25 * rho.inj^0.75
        dI.prod <- 0.241 / Dw^4.75 * mu.prod^0.25 * rho.prod^0.75
        I.well <- I.w0 * (a.inj / N.inj^1.75 + 1 / N.prod^1.75) + (z[i]-Z0) * (dI.inj / N.inj + dI.prod / N.prod)
        eta.el <- 0.078795*log(h.prod)-1.00081      #Zarrouk & Moon, 2012
        #计算V.opt
        dV.opt <- calc_dV.opt(z[i], I.well, I.res, N.res, eta.el, rho.inj, rho.prod, cw, T0, T.gradient, T.inj, g)
        #分配得到V.cap
        dV.max <- calc_dV.max(dV.opt, sigma, z[i], I.res, g, T.inj)
        
        E.th <- rho.inj * dV.max * (h.prod-h.inj)
        DP.g <- z[i]*g * (rho.inj-rho.prod)
        W.EGS <- (I.well * dV.max^1.75 + I.res/N.res * dV.max - DP.g) * dV.max
        E.el[i,rdm] <- eta.el * E.th - W.EGS
        #4:9这个深度各算1万遍
        #到此为止，产生的电能计算结束




        #开始计算价格
        #钻井最小价格
        C.well.min <- (1.72e-7*z[i]^2+2.3e-3*z[i]-0.62)*1e6	  #USD, z in m
        #钻井价格 4:9, 10000次
        C.well[i,rdm] <- runif(1, min=C.well.min, max=2*C.well.min)       #USD/well
        #电厂建设最小价格
        C.plant.min <- 750+1125*exp(-0.006115*(E.el[i,rdm]*1e-6 - 5))		#Pnet in MW, cost USD/kW p. 7-14 (251) MIT 2006
        #电厂建设价格，4:9,10000次
        C.plant <- runif(1, min=C.plant.min, max=2*C.plant.min)
        
        Dt.EGS[i,rdm] <- runif(1, min=Dt.EGS.min, max=Dt.EGS.max)
        Dt.well <- runif(1, min=Dt.well.min, max=Dt.well.max)
        #round四舍五入函数
        N.well.renew <- round((Dt.EGS[i,rdm]-Dt.well)/Dt.well)
        indneg <- which(N.well.renew < 0)
        N.well.renew[indneg] <- 0
        #基本价格
        C.EGS[i,rdm] <- N.well*(C.well[i,rdm]+C.frac)+C.plant*E.el[i,rdm]*1e-3 +	N.well*N.well.renew*C.well[i,rdm]	#USD
        E.EGS[i,rdm] <- E.el[i,rdm]*1e-3*n.hr*Dt.EGS[i,rdm]
        P.EGS[i,rdm] <- C.EGS[i,rdm]/E.EGS[i,rdm]
        #基本的价格在这里计算完成

        #下面开始计算heat credit
        #heat credit
        #heat credit 和L有关，所以需要进行三维
        E.DH0 <- (1-eta.el)*E.th
        eta.DH <- exp(W0(-U*L*DTs/E.DH0))    #nL vector
        ind0 <- which(E.DH0 < exp(1)*U*L*DTs)
        eta.DH[ind0] <- NA
        dV.DHnom <- -U*L/(rho.T(Ts)*cw*W0(-U*L*DTs/E.DH0))
        ind00 <- which(dV.DHnom < dV.max)
        eta.DH[ind00] <- NA
        #大小是 6*10000*81
        E.DH[i,rdm,] <- eta.DH * E.DH0
        
        C.DH[i,rdm,] <- (C1+C2*D.DH)*L  *1.2    #EUR to USD
        En.DH[i,rdm,] <- E.DH[i,rdm,]*1e-3*2500*Dt.EGS[i,rdm]
        E.comb[i,rdm,] <- E.EGS[i,rdm] + En.DH[i,rdm,]/3
        C.comb[i,rdm,] <- C.EGS[i,rdm] + C.DH[i,rdm,]
        P.comb[i,rdm,] <- C.comb[i,rdm,]/E.comb[i,rdm,]
    }
}


#在这里得到价格
# nz <- 6, nL <- 10000
Pnet_z <- rep(NA, nz)
costs.well_z <- rep(NA, nz)
costs_z <- rep(NA, nz)
price_z <- rep(NA, nz)

E.combL0_z <- array(NA, dim=c(nz,nL))
costs.combL0_z <- array(NA, dim=c(nz,nL))
price.combL0_z <- array(NA, dim=c(nz,nL))

for(i in 1:nz) {
  Pnet_z[i] <- median(E.el[i, ])
  costs.well_z[i] <- median(C.well[i, ])
  costs_z[i] <- median(C.EGS[i, ])
  price_z[i] <- median(P.EGS[i, ])
  
  for(j in 1:nL) {
    E.combL0_z[i,j] <- median(E.comb[i, ,j], na.rm=T)
    costs.combL0_z[i,j] <- median(C.comb[i, ,j], na.rm=T)
    price.combL0_z[i,j] <- median(P.comb[i, ,j], na.rm=T)
  }
}

#至此，价格全部结束。没有加入死亡的概率的价格是不完美的价格