library(ggplot2)
rho.T <- function( T ) return(1030 - 0.1625 * T - 0.00269 * T ^ 2 )
h.T <- function(T, cw) return((T+273)*cw)   
mu.T <- function(T) return(343.18e-7*10^(247.8/(T.inj+140))) 
T0 <- 15 #deg C
cw <- 4180
I.res <- 0.2*1e9
a.inj <- 2
Z0 <- 2500           #m
Dw <- 0.25           #m
g <- 9.81            #m/s2

type <- "triplet"
z <- seq(4,9,1)*1e3
I.w0 <- 1869*(1e3)^1.75     #Pa/(m3/s)^1.75   to check (Basel)
T.gradient.min <- 30       #deg C
T.gradient.max <- 40       #deg C
T.gradient <- c(3, 40)
T.inj.min <- 60                   #deg C
T.inj.max <- 75                   #deg C
dT.dZ <- c(30,40)
V.inj <- seq(0, 0.4, 0.001)

T.inj.min <- 60                   #deg C
T.inj.max <- 75                   #deg C

if(type == "doublet") {N.res <- 1; N.inj <- 1; N.prod <- 1; N.well <- 2}
if(type == "triplet") {N.res <- 2; N.inj <- 2; N.prod <- 1; N.well <- 3}
if(type == "3stageDoublet") {N.res <- 3; N.inj <- 1; N.prod <- 1; N.well <- 2}
E.el <-array(NA, dim = c(6,length(V.inj),2))
for(i in 1:2){
  for(j in 1:6){
    for( rdm in 1:length(V.inj)){
      
      T.inj <- runif(1, min=T.inj.min, max=T.inj.max)
      
      T.prod <- T0+T.gradient[i]*z[j]*1e-3
      h.inj <- h.T(T.inj, cw)
      h.prod <- h.T(T.prod, cw)
      rho.inj <- rho.T(T.inj)
      rho.prod <- rho.T(T.prod)
      mu.inj <- mu.T(T.inj)
      mu.prod <- mu.T(T.prod)
      
      DT <- T0+ z[j]/1000*dT.dZ[i] - T.inj
      
      E.th <- rho.inj*V.inj[rdm]*cw*DT
      ## DP.res ok
      DP.res <- I.res*V.inj[rdm]/N.res
      ## DP.well ok
      dI.inj <- 0.241 / Dw^4.75 * mu.inj^0.25 * rho.inj^0.75
      dI.prod <- 0.241 / Dw^4.75 * mu.prod^0.25 * rho.prod^0.75
      I.well <- I.w0 * (a.inj / N.inj^1.75 + 1 / N.prod^1.75) + (z[i]-Z0) * (dI.inj / N.inj + dI.prod / N.prod)
      DP.well <- I.well*V.inj[rdm]^1.75
      
      DP.g <- z[j]*g*(rho.inj-rho.prod)
      W.EGS <-(DP.res + DP.well - DP.g) * V.inj[rdm]
      ##eta.el ok
      eta.el <- 0.078795*log(h.prod)-1.00081      #Zarrouk & Moon, 2012
      
      E.el[j,rdm,i] = eta.el*E.th - W.EGS
    }
  }
}


##在这里开始绘图

type <- rep(c(4:9),each = length(V.inj))
x_value <- rep(seq(0,0.4,0.001),times = 6)
y <- expression(dot(E[el]))
y_value1 <- c(E.el[1,,1],E.el[2,,1],E.el[3,,1],E.el[4,,1],E.el[5,,1],E.el[6,,1])
E1 <- data.frame(x = x_value, y = y_value1,type = type)
#here I need to modify some parameters to make the figure beautiful
new_E <- E1[y_value1>0,]
g1 <- ggplot(data = new_E, mapping = aes(x = x, y = y/1e6,group = type ,type = type ,color = type))+
geom_smooth()+
theme(plot.title = element_text(hjust = 0.5))+
labs(title = "Triplet,dT/dz=30K/km", x = expression(dot(V)(m^3/s)), y = y, lty = "Split")+
scale_y_continuous(breaks = c(2,4,6,8,10,20,40,60))


y_value2 <- c(E.el[1,,2],E.el[2,,2],E.el[3,,2],E.el[4,,2],E.el[5,,2],E.el[6,,2])
E2 <- data.frame(x = x_value, y = y_value2,type = type)
#here I need to modify some parameters to make the figure beautiful
new_E <- E2[y_value2>0,]

g2 <- ggplot(data = new_E, mapping = aes(x = x, y = y/1e6,group = type ,type = type ,color = type))+
geom_smooth()+
theme(plot.title = element_text(hjust = 0.5))+
labs(title = "Triplet,dT/dz=40K/km", x = expression(dot(V)(m^3/s)), y = y,lty = "Split")+
scale_y_continuous(breaks = c(2,4,6,8,10,20,40,60))

grid.arrange(g1,g2,ncol=1,nrow=2)