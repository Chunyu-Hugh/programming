#FUNCTIONS
library(ggplot2)



rho.T <- function(T) return(1030-0.1625*T-0.00269*T^2) #kg/m3
h.T <- function(T, cw) return((T+273)*cw)   #J/kg, T in degree Celsius
mu.T <- function(T) return(243.18e-7*10^(247.8/(T.inj+273-140))) #viscosity McDermott et al.(2006)
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


#input parameters EGS
z <- seq(4, 9, .1)*1e3         #depth m
type <- "triplet"
T.gradient <- 40       #deg C
T.inj <- 60                   #deg C
T0 <- 15                      #deg C
cw <- 4180           #J/(kg.K)
g <- 9.81            #m/s2
I.res <- 1e8           #Pa/(m3/s) to check
I.w0 <- 1869*(1e3)^1.75     #Pa/(m3/s)^1.75   to check (Basel)
a.inj <- 2
Z0 <- 2500           #m
Dw <- 0.25           #m
sigma <- 0.1

#input district heating
L <- seq(0,300,1)*1e3    #m
U <- 8                     #W/m/K
Ts <- 50                 #supply temp
DTs <- Ts-T0


nz <- length(z)
nL <- length(L)

E.el <- rep(NA, nz)
dV.opt <-  rep(NA, nz)
dV.max <-  rep(NA, nz)
eta.el <-  rep(NA, nz)
W.EGS <- rep(NA, nz)
E.th <- rep(NA, nz)
E.el <-  rep(NA, nz)
for(i in 1:nz) {

  if(type == "doublet") {N.res <- 1; N.inj <- 1; N.prod <- 1; N.well <- 2}
  if(type == "triplet") {N.res <- 2; N.inj <- 2; N.prod <- 1; N.well <- 3}
  if(type == "3stageDoublet") {N.res <- 3; N.inj <- 1; N.prod <- 1; N.well <- 2}

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
  eta.el[i] <- 0.078795*log(h.prod)-1.00081      #Zarrouk & Moon, 2012
  
  dV.opt[i] <- calc_dV.opt(z[i], I.well, I.res, N.res, eta.el[i], rho.inj, rho.prod, cw, T0, T.gradient, T.inj, g)
  dV.max[i] <- calc_dV.max(dV.opt[i], sigma, z[i], I.res, g, T.prod)
  
  E.th[i] <- rho.inj * dV.max[i] * (h.prod-h.inj)
  DP.g <- z[i]*g * (rho.inj-rho.prod)
  W.EGS[i] <- (I.well * dV.max[i]^1.75 + I.res/N.res * dV.max[i] - DP.g) * dV.max[i]
  E.el[i] <- eta.el[i] * E.th[i] - W.EGS[i]
}
E.DH0 <- (1-eta.el)*E.th
#此时就可以画图啦，画E.th,E.DH0,E.el
dghy <- data.frame(x=rep(z,3)*1e-3, y=c(E.th,E.el,E.DH0)*1e-6, id=c(rep("Eth",nz), rep("Eel",nz), rep("EDH,0",nz)))
ggplot(data = dghy)+
geom_line(aes(x=x,y=y,group = id, lty = id))+
theme_minimal()+
labs(title = "(a)", x = "z (km)", y = "Power (MW)", lty="Split")

#在这里计算eta.DH
eta.DH <- array(NA, dim = c(nz, nL))
eta.DH.eco <- array(NA, dim = c(nz, nL))
E.DH <- array(NA, dim = c(nz, nL))
E.DH.eco <- array(NA, dim = c(nz, nL))

for(i in 1:nz) {
  eta.DH[i,] <- exp(W0(-U*L*DTs/E.DH0[i]))
  ind0 <- which(E.DH0[i] < exp(1)*U*L*DTs)
  eta.DH.eco[i,] <- eta.DH[i,]
  eta.DH.eco[i, ind0] <- NA
  
  dV.DHnom <- -U*L/(rho.T(Ts)*cw*W0(-U*L*DTs/E.DH0[i]))
  ind00 <- which(dV.DHnom < dV.max[i])
  eta.DH.eco[i, ind00] <- NA
  
  E.DH[i,] <- eta.DH[i,] * E.DH0[i]
  E.DH.eco[i,] <- eta.DH.eco[i,] * E.DH0[i]
}

grid <- expand.grid(x=L, y=z)
grid.ind <- expand.grid(i=seq(nL), j=seq(nz))
npt <- nrow(grid) #15351个数字
E.DH.gg <- rep(NA, npt)
E.DH.eco.gg <- rep(NA, npt)
E.el.gg <- rep(NA, npt)
for(pt in 1:npt) {
  j <- grid.ind$i[pt]
  i <- grid.ind$j[pt]
  E.DH.gg[pt] <- E.DH[i,j]
  E.DH.eco.gg[pt] <- E.DH.eco[i,j]
  E.el.gg[pt] <- E.el[i]
}
DH.gg <- data.frame(grid, E.DH=E.DH.gg*1e-6, E.DH.eco=E.DH.eco.gg*1e-6, E.el=E.el.gg*1e-6)




##下面是画图
pdf(file="fig3.pdf",family="Times",pointsize=16)
#get!
g1 <- ggplot(data=data.frame(x=rep(z,3)*1e-3, y=c(E.th,E.el,E.DH0)*1e-6, id=c(rep("Eth",nz), rep("Eel",nz), rep("EDH,0",nz)))) +
  geom_line(aes(x=x, y=y, group=id, lty=id)) +
  theme_minimal() +
  labs(title = "(a)", x = "z (km)", y = "E (MW)", lty="Split")
#get!
zi <- seq(4,9)
nzi <- length(zi)
g2 <- ggplot() +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH[which(z == zi[1]*1e3),], z=rep(zi[1],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH[which(z == zi[2]*1e3),], z=rep(zi[2],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH[which(z == zi[3]*1e3),], z=rep(zi[3],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH[which(z == zi[4]*1e3),], z=rep(zi[4],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH[which(z == zi[5]*1e3),], z=rep(zi[5],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH[which(z == zi[6]*1e3),], z=rep(zi[6],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH.eco[which(z == zi[1]*1e3),], z=rep(zi[1],nL)), aes(x=x, y=y, col=z)) +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH.eco[which(z == zi[2]*1e3),], z=rep(zi[2],nL)), aes(x=x, y=y, col=z)) +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH.eco[which(z == zi[3]*1e3),], z=rep(zi[3],nL)), aes(x=x, y=y, col=z)) +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH.eco[which(z == zi[4]*1e3),], z=rep(zi[4],nL)), aes(x=x, y=y, col=z)) +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH.eco[which(z == zi[5]*1e3),], z=rep(zi[5],nL)), aes(x=x, y=y, col=z)) +
  geom_line(data=data.frame(x=L*1e-3, y=eta.DH.eco[which(z == zi[6]*1e3),], z=rep(zi[6],nL)), aes(x=x, y=y, col=z)) +
  scale_colour_gradient(low = "#64641c", high = "red") +
  theme_minimal() +
  labs(title = "(b)", x = "L (km)", y = expression(eta[DH]), col="z (km)")
#get!
g3 <- ggplot() +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH[which(z == zi[1]*1e3),]*1e-6, z=rep(zi[1],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH[which(z == zi[2]*1e3),]*1e-6, z=rep(zi[2],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH[which(z == zi[3]*1e3),]*1e-6, z=rep(zi[3],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH[which(z == zi[4]*1e3),]*1e-6, z=rep(zi[4],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH[which(z == zi[5]*1e3),]*1e-6, z=rep(zi[5],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH[which(z == zi[6]*1e3),]*1e-6, z=rep(zi[6],nL)), aes(x=x, y=y, col=z), lty="dotted") +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH.eco[which(z == zi[1]*1e3),]*1e-6, z=rep(zi[1],nL)), aes(x=x, y=y, col=z)) +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH.eco[which(z == zi[2]*1e3),]*1e-6, z=rep(zi[2],nL)), aes(x=x, y=y, col=z)) +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH.eco[which(z == zi[3]*1e3),]*1e-6, z=rep(zi[3],nL)), aes(x=x, y=y, col=z)) +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH.eco[which(z == zi[4]*1e3),]*1e-6, z=rep(zi[4],nL)), aes(x=x, y=y, col=z)) +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH.eco[which(z == zi[5]*1e3),]*1e-6, z=rep(zi[5],nL)), aes(x=x, y=y, col=z)) +
  geom_line(data=data.frame(x=L*1e-3, y=E.DH.eco[which(z == zi[6]*1e3),]*1e-6, z=rep(zi[6],nL)), aes(x=x, y=y, col=z)) +
  scale_colour_gradient(low = "yellow", high = "red") +
  theme_minimal() +
  labs(title = "(c)", x = "L (km)", y = expression(paste(E[DH], " (MW)"), sep=""), col="z (km)")

g4 <- ggplot(data=DH.gg, aes(x=x*1e-3, y=-y*1e-3)) +
  geom_raster(aes(fill=E.DH.eco)) +
  scale_fill_gradient(low="yellow", high="red", limits=c(0,255), na.value="lightgrey") +
  geom_contour(aes(z=E.DH), breaks = 10^seq(1,2.3,0.1), col="black", lty="dotted") +
  theme_minimal() +
#  theme(legend.position="none") +
  labs(title = "(d)", x = "L (km)", y = "z (km)", fill=expression(paste(E[DH], " (MW)"), sep=""))

grid.arrange(g1,g2,g3, ncol=1, nrow=3)
dev.off()





print("DONE")