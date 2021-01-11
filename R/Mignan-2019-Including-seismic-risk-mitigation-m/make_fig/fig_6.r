library(ggplot2)
x <- c(-2.6, -3.2, -2.0, -1.4, -2.4, -3.8, - 3.1, -0.5, -0.9, 0.1, -4.2, -2.8, -1.6)
y <- c(0.7,  0.8,  1.4, 0.9, 1.1, 2.2, 1.8, 1.1, 0.8, 1.6, 1.1, 0.8 ,1.0)
ggplot(data = data.frame(a = x, b = y), aes(x=a, y=b))+
geom_point()+
labs(title = "(a) Underground feedback", x = expression(a(m^-3)))+
theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行




#b图
library(ggplot2)

x2 <- seq(0,1,0.001)*1e8
y2 <- 2*log10(3.10e10*x2)/3 -10.7 +14/3
ggplot(data = data.frame(x2 = x2,y2 = y2),aes(x=x2,y=y2))+
geom_line()


#c图

#California
IPE1 <- function(m,r){
  #Atkinson and Wald, 2007 - California
  c1 <- 12.27
  c2 <- 2.270
  c3 <- 0.1304
  c4 <- -1.30
  c5 <- -0.0007070
  c6 <- 1.95
  c7 <- -0.577
  Rt <- 30
  sigma <- 0.4
  int <- c1+c2*(m-6)+c3*(m-6)^2+c4*log10(r)+c5*r+c6*B+c7*m*log10(r)
  return(list(median=int, sigma=sigma))
}
#Central/Eastern U.S.A
IPE2 <- function(m,r){
  #Atkinson and Wald, 2007 - Central/Eastern US
  c1 <- 11.72
  c2 <- 2.36
  c3 <- 0.1155
  c4 <- -0.44
  c5 <- -0.002044
  c6 <- 2.31
  c7 <- -0.479
  Rt <- 80
  sigma <- 0.4
  int <- c1+c2*(m-6)+c3*(m-6)^2+c4*log10(r)+c5*r+c6*B+c7*m*log10(r)
  return(list(median=int, sigma=sigma))
}
#Switzerland
IPE3 <- function(m,d){
  # m moment magnitude, d hypocentral distance (km)
  # int observed Mercalli Intensity
  a <- -0.69182; b <- 0.00084; alpha <- 0.7317; beta <- 1.2567
  c0 <- beta; c1 <- alpha; c2 <- -alpha*a; c3 <- -alpha*b
  int <- (m-c2*log(d/30)-c3*(d-30)-c0)/c1
  return(list(median=, sigma=1))
} 
#Global
IPE4 <- function(m,d){
  # Allen et al. (2012), pp. 418, Eqs. 2 and 3 (here without site effect), pp. 423, Eq. 6
  # m moment magnitude, d hypocentral distance (km)
  # int generic Cancani scale (12 degrees)
  c0 <- 2.085; c1 <- 1.428; c2 <- -1.402; c4 <- 0.078
  m1 <- -0.209; m2 <- 2.042
  Rm <- m1+m2*exp(m-5)
  
  #	if(d <= 50) int <- c0+c1*m+c2*log(sqrt(d^2+Rm^2))
  #	else int <- c0+c1*m+c2*log(sqrt(d^2+Rm^2))+c4*log(d/50)
  
  int <- c0+c1*m+c2*log(sqrt(d^2+Rm^2))+c4*log(d/50)
  indnear <- which(d <= 50)
  if(length(indnear) != 0) int[indnear] <- c0+c1*m+c2*log(sqrt(d[indnear]^2+Rm^2))
  
  s1 <- 0.82; s2 <- 0.37; s3 <- 22.9
  sigma <- s1+s2/(1+(d/s3)^2)
  return(list(median=int, sigma=sigma))
}

for(i in 1:)