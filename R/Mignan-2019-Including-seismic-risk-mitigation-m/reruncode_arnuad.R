# 文章的主要目的是计算LCOE、计算最优选址�?
# rerun the code
rm(list = ls())
#LIBRARIES
library(ggplot2)
library(gridExtra) #grid.arrange()
library(reshape2) #melt()
library(directlabels)
library(msm) #ptnorm()
library(signal) #interp1()
library(RColorBrewer)
library(splancs) #inpip()


#SETUP
#使其在R这个文件夹下面进行运�?
wd <- paste(getwd(), sep = "")
wdsetwd(wd)
outd <- "outputs"
if (!file.exists(outd)) dir.create(outd)
figd <- "figures"
if (!file.exists(figd)) dir.create(figd)

# 经济部分的计�?

##################################### ECONOMIC PART #########################################

#FUNCTIONS
rho.T <- function(T) return(1030 - 0.1625 * T - 0.00269 * T ^ 2) #kg/m3
h.T <- function(T, cw) return((T + 273) * cw) #J/kg, T in degree Celsius摄氏温度
mu.T <- function(T) return(343.18e-7 * 10 ^ (247.8 / (T.inj + 140))) #viscosity McDermott et al. (2006)
# 计算最值，最合适的�?

calc_dV.opt <- function(z, I.well, I.res, N.res, eta.el, rho.inj, rho.prod, cw, T0, T.gradient, T.inj, g) {
  a <- -I.well
  b <- -I.res / N.res
  c <- eta.el * rho.inj * cw * (T0 + z * 1e-3 * T.gradient - T.inj) + z * g * (rho.inj - rho.prod)

  a1 <- 2.75 * a #  a1 <- 1/4*2.75*a
  b1 <- 2 * b
  c1 <- c
  V0.p <- (-b1 + sqrt(b1 ^ 2 - 4 * a1 * c1)) / (2 * a1)
  V0.m <- (-b1 - sqrt(b1 ^ 2 - 4 * a1 * c1)) / (2 * a1)

  if (V0.p > 0) V0 <- V0.p else V0 <- V0.m
  a2 <- 1 / 2 * 2.75 * 1.75 * 0.75 * a * V0 ^ (-0.25)
  b2 <- 2.75 * 1.75 * a * V0 ^ (0.75) + 2 * b
  c2 <- 2.75 * a * V0 ^ (1.75) + 2 * b * V0 + c
  DV.p <- (-b2 + sqrt(b2 ^ 2 - 4 * a2 * c2)) / (2 * a2)
  DV.m <- (-b2 - sqrt(b2 ^ 2 - 4 * a2 * c2)) / (2 * a2)
  if (DV.p > 0) DV <- DV.p else DV <- DV.m

  Vopt <- V0 + DV
  return(Vopt)
}

# 这是有一个公式，分段函数寻找需要哪一个值进行计�?
calc_dV.max <- function(dV.opt, sigma, z, I.res, g, Temp) {
  nu <- 0.25
  Sv <- 2500 * g * z
  P.hydro <- rho.T(Temp) * g * z
  ratio <- nu / (1 - nu)
  #Smin-P.hydro = nu/(1-nu) (Sv-P.hydro)
  Smin <- (1 - ratio) * P.hydro + ratio * Sv

  dV.max <- (Smin * (1 - 3 * sigma) - P.hydro) / I.res
  if (dV.max > dV.opt) dV.max <- dV.opt
  return(dV.max)
}

W0 <- function(y) {
  x <- y - y ^ 2 + 3 / 2 * y ^ 3 - 8 / 3 * y ^ 4 + 125 / 24 * y ^ 5
  return(x)
}
# 出来的是圆的数据
circle <- function(center, radius, npts = 100) {
  pts <- seq(0, 2 * pi, length.out = npts)
  x <- center[1] + radius * cos(pts)
  y <- center[2] + radius * sin(pts)
  return(data.frame(x = x, y = y))
}
# 输入参数EGS
#input parameters EGS
T0 <- 15 #deg C
cw <- 4180 #J/(kg.K)
g <- 9.81 #m/s2
I.res <- 0.2 * 1e9

I.w0 <- 1869 * (1e3) ^ 1.75 #Pa/(m3/s)^1.75   to check (Basel)
a.inj <- 2
Z0 <- 2500 #m
Dw <- 0.25 #m
sigma <- 0 #0.05   #0
U <- 8 #W/m/K
Ts <- 50 #supply temp
DTs <- Ts - T0
n.hr <- 8000 #hrs/yr
C.frac <- 1e6 #USD/well
C1 <- 151 #EUR/m
C2 <- 1378 #EUR/m2
D.DH <- 0.5 #m

type <- "triplet"

z <- seq(4, 9, 1) * 1e3 #depth m
L <- seq(0, 80, 1) * 1e3 #m

T.gradient.min <- 30 #deg C
T.gradient.max <- 40 #deg C
T.inj.min <- 60 #deg C
T.inj.max <- 75 #deg C
Dt.EGS.min <- 20 #yr
Dt.EGS.max <- 30 #yr
Dt.well.min <- 5 #yr
Dt.well.max <- 30 #yr


nz <- length(z)
nL <- length(L)
nrdm <- 1e4
# nz <- 6
# nL <- 81
# nrdm <- 10000
# 6�?10000�?
E.el <- array(NA, dim = c(nz, nrdm))
Dt.EGS <- array(NA, dim = c(nz, nrdm))
C.well <- array(NA, dim = c(nz, nrdm))
C.EGS <- array(NA, dim = c(nz, nrdm))
E.EGS <- array(NA, dim = c(nz, nrdm))
P.EGS <- array(NA, dim = c(nz, nrdm))
# 81�? 6�?10000�?
E.DH <- array(NA, dim = c(nz, nrdm, nL))
C.DH <- array(NA, dim = c(nz, nrdm, nL))
En.DH <- array(NA, dim = c(nz, nrdm, nL))
E.comb <- array(NA, dim = c(nz, nrdm, nL))
C.comb <- array(NA, dim = c(nz, nrdm, nL))
P.comb <- array(NA, dim = c(nz, nrdm, nL))

for (i in 1:nz) {
  for (rdm in 1:nrdm) {
    # 每一个距离都�?1000个数 
    # 1000�?
    # 直接选取triplet
    if (type == "doublet") { N.res <- 1; N.inj <- 1; N.prod <- 1; N.well <- 2 }
    if (type == "triplet") { N.res <- 2; N.inj <- 2; N.prod <- 1; N.well <- 3 }
    if (type == "3stageDoublet") { N.res <- 3; N.inj <- 1; N.prod <- 1; N.well <- 2 }
    # 均匀分布
    T.gradient <- runif(1, min = T.gradient.min, max = T.gradient.max)
    T.inj <- runif(1, min = T.inj.min, max = T.inj.max)

    T.prod <- T0 + T.gradient * z[i] * 1e-3
    h.inj <- h.T(T.inj, cw)
    h.prod <- h.T(T.prod, cw)
    rho.inj <- rho.T(T.inj)
    rho.prod <- rho.T(T.prod)
    mu.inj <- mu.T(T.inj)
    mu.prod <- mu.T(T.prod)

    dI.inj <- 0.241 / Dw ^ 4.75 * mu.inj ^ 0.25 * rho.inj ^ 0.75
    dI.prod <- 0.241 / Dw ^ 4.75 * mu.prod ^ 0.25 * rho.prod ^ 0.75
    I.well <- I.w0 * (a.inj / N.inj ^ 1.75 + 1 / N.prod ^ 1.75) + (z[i] - Z0) * (dI.inj / N.inj + dI.prod / N.prod)
    eta.el <- 0.078795 * log(h.prod) - 1.00081 #Zarrouk & Moon, 2012
    #计算 V.opt
    dV.opt <- calc_dV.opt(z[i], I.well, I.res, N.res, eta.el, rho.inj, rho.prod, cw, T0, T.gradient, T.inj, g)
    #分配V.cap
    dV.max <- calc_dV.max(dV.opt, sigma, z[i], I.res, g, T.inj)

    E.th <- rho.inj * dV.max * (h.prod - h.inj)
    DP.g <- z[i] * g * (rho.inj - rho.prod)
    W.EGS <- (I.well * dV.max ^ 1.75 + I.res / N.res * dV.max - DP.g) * dV.max

    #第i�? 第rdm列的值是�?
    E.el[i, rdm] <- eta.el * E.th - W.EGS


    # 计算出来的电能也满足均匀分布



    #开始计算价�?
    #钻井价格的数�? 第i�? 第rdm�? 是按照均匀分布的随机数
    C.well.min <- (1.72e-7 * z[i] ^ 2 + 2.3e-3 * z[i] - 0.62) * 1e6 #USD, z in m
    C.well[i, rdm] <- runif(1, min = C.well.min, max = 2 * C.well.min) #USD/well

    #建厂价格是按照均匀分布的随机数
    C.plant.min <- 750 + 1125 * exp(-0.006115 * (E.el[i, rdm] * 1e-6 - 5)) #Pnet in MW, cost USD/kW p. 7-14 (251) MIT 2006
    C.plant <- runif(1, min = C.plant.min, max = 2 * C.plant.min)

    Dt.EGS[i, rdm] <- runif(1, min = Dt.EGS.min, max = Dt.EGS.max)

    Dt.well <- runif(1, min = Dt.well.min, max = Dt.well.max)

    N.well.renew <- round((Dt.EGS[i, rdm] - Dt.well) / Dt.well)
    indneg <- which(N.well.renew < 0)
    N.well.renew[indneg] <- 0

    # 未加入其余风险的价格，也是服从均匀分布的数字啦
    C.EGS[i, rdm] <- N.well * (C.well[i, rdm] + C.frac) + C.plant * E.el[i, rdm] * 1e-3 + N.well * N.well.renew * C.well[i, rdm] #USD
    E.EGS[i, rdm] <- E.el[i, rdm] * 1e-3 * n.hr * Dt.EGS[i, rdm]
    P.EGS[i, rdm] <- C.EGS[i, rdm] / E.EGS[i, rdm]



    # 计算产热 heat credit
    # 产热不是均匀分布的随机数，是有序的数�?
    E.DH0 <- (1 - eta.el) * E.th
    eta.DH <- exp(W0(-U * L * DTs / E.DH0)) #L 向量
    ind0 <- which(E.DH0 < exp(1) * U * L * DTs) # 其实可以使用判断，但是没�?
    eta.DH[ind0] <- NA # ind0 是返回的坐标

    dV.DHnom <- -U * L / (rho.T(Ts) * cw * W0(-U * L * DTs / E.DH0))
    ind00 <- which(dV.DHnom < dV.max)
    eta.DH[ind00] <- NA # ind0是返回的坐标
    # 81�? 6�?10000�?
    E.DH[i, rdm,] <- eta.DH * E.DH0
    C.DH[i, rdm,] <- (C1 + C2 * D.DH) * L * 1.2 #EUR to USD
    En.DH[i, rdm,] <- E.DH[i, rdm,] * 1e-3 * 2500 * Dt.EGS[i, rdm]

    E.comb[i, rdm,] <- E.EGS[i, rdm] + En.DH[i, rdm,] / 3
    C.comb[i, rdm,] <- C.EGS[i, rdm] + C.DH[i, rdm,]
    # 价格合并
    # 价格也是随机数，需要进行抽�?
    P.comb[i, rdm,] <- C.comb[i, rdm,] / E.comb[i, rdm,]
  }
}
# nz <- 6
# 6�? 深度�?4000 5000 6000 7000 8000 9000

######## 完全看懂的时�?########################
######## 2021�?2�?25�?15�?41�?#################



############# next part ######################
# nz <- 6
# nL <- 81
# 这四个值每一个进行复制，复制六次
Pnet_z <- rep(NA, nz)
costs.well_z <- rep(NA, nz)
costs_z <- rep(NA, nz)
price_z <- rep(NA, nz)
# 6*81的数�?
E.combL0_z <- array(NA, dim = c(nz, nL))
costs.combL0_z <- array(NA, dim = c(nz, nL))
price.combL0_z <- array(NA, dim = c(nz, nL))
for (i in 1:nz) {
  # median() 计算下面这些数字每一行的中位�?
  Pnet_z[i] <- median(E.el[i,])
  costs.well_z[i] <- median(C.well[i,])
  costs_z[i] <- median(C.EGS[i,])
  price_z[i] <- median(P.EGS[i,])
  for (j in 1:nL) {
    # 因为�?81�? 所以每一个的 每一个距离的行都需要取中位�?
    # 6�?81�?
    E.combL0_z[i, j] <- median(E.comb[i,, j], na.rm = T)
    costs.combL0_z[i, j] <- median(C.comb[i,, j], na.rm = T)
    price.combL0_z[i, j] <- median(P.comb[i,, j], na.rm = T)
  }
}

#################### 计算LCOE ####################
z <- seq(4, 9, 1)
#L转换为km为单位制�?
d <- L * 1e-3 #km
nd <- length(d) # 81 管道
nz <- length(z) # 6  钻井距离的区�?
#grid 的格式是 第一列是管道距离 第二列是打井深度
grid <- expand.grid(x = d, y = z) # 486�? 两列

grid.ind <- expand.grid(i = seq(nd), j = seq(nz))
npt <- nrow(grid)
# npt <- 486
LCOE0.gg <- rep(NA, npt)
LCOEcredit.gg <- rep(NA, npt)
# npt <- 486
for (pt in 1:npt) {
  j <- grid.ind$i[pt]
  i <- grid.ind$j[pt]
  LCOE0.gg[pt] <- price_z[i]
  LCOEcredit.gg[pt] <- price.combL0_z[i, j]
}
# 合并 得到价格
LCOE.gg <- data.frame(grid, P_standard = LCOE0.gg * 1e2, P_heatcredit = LCOEcredit.gg * 1e2)



##################################### RISK PART #########################################

#FUNCTIONS
IPE1 <- function(m, r) {
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
  B <- rep(0, length(r))
  indfar <- which(r > Rt)
  if (length(indfar) != 0) B[indfar] <- log10(r[indfar] / Rt)
  int <- c1 + c2 * (m - 6) + c3 * (m - 6) ^ 2 + c4 * log10(r) + c5 * r + c6 * B + c7 * m * log10(r)
  return(list(median = int, sigma = sigma))
}

IPE2 <- function(m, r) {
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
  B <- rep(0, length(r))
  indfar <- which(r > Rt)
  if (length(indfar) != 0) B[indfar] <- log10(r[indfar] / Rt)
  int <- c1 + c2 * (m - 6) + c3 * (m - 6) ^ 2 + c4 * log10(r) + c5 * r + c6 * B + c7 * m * log10(r)
  return(list(median = int, sigma = sigma))
}

IPE3 <- function(m, d) {
  # m moment magnitude, d hypocentral distance (km)
  # int observed Mercalli Intensity
  a <- -0.69182;
  b <- 0.00084;
  alpha <- 0.7317;
  beta <- 1.2567
  c0 <- beta;
  c1 <- alpha;
  c2 <- -alpha * a;
  c3 <- -alpha * b
  return(list(median = (m - c2 * log(d / 30) - c3 * (d - 30) - c0) / c1, sigma = 1))
}

IPE4 <- function(m, d) {
  # Allen et al. (2012), pp. 418, Eqs. 2 and 3 (here without site effect), pp. 423, Eq. 6
  # m moment magnitude, d hypocentral distance (km)
  # int generic Cancani scale (12 degrees)
  c0 <- 2.085;
  c1 <- 1.428;
  c2 <- -1.402;
  c4 <- 0.078;
  m1 <- -0.209;
  m2 <- 2.042
  Rm <- m1 + m2 * exp(m - 5)

  #	if(d <= 50) int <- c0+c1*m+c2*log(sqrt(d^2+Rm^2))
  #	else int <- c0+c1*m+c2*log(sqrt(d^2+Rm^2))+c4*log(d/50)

  int <- c0 + c1 * m + c2 * log(sqrt(d ^ 2 + Rm ^ 2)) + c4 * log(d / 50)
  indnear <- which(d <= 50)
  if (length(indnear) != 0) int[indnear] <- c0 + c1 * m + c2 * log(sqrt(d[indnear] ^ 2 + Rm ^ 2))

  s1 <- 0.82;
  s2 <- 0.37;
  s3 <- 22.9
  sigma <- s1 + s2 / (1 + (d / s3) ^ 2)
  return(list(median = int, sigma = sigma))
}
# 风险曲线
fct_hazardcurve <- function(rate, intensity, sigma, nsigma, increment, distrib) {
  int <- increment;
  nbi <- length(int)
  nbm <- length(rate)
  rate_mi <- numeric(nbm * nbi);
  dim(rate_mi) <- c(nbm, nbi)
  rate_all <- numeric(nbi)
  for (i in 1:nbm) {
    if (distrib == "deterministic") {
      Pr_exceed <- rep(0, nbi)
      Pr_exceed[which(int < intensity[i])] <- 1
    }
    #		if(distrib == "lognormal") Pr_exceed <- 1-plnorm(int, meanlog=log(intensity[i]), sdlog=sigma[i])
    #		if(distrib == "normal") Pr_exceed <- 1-pnorm(int, mean=intensity[i], sd=sigma[i])
    #		indlow <- which(int < intensity[i]-3*sigma[i])
    #		indhigh <- which(int > intensity[i]+3*sigma[i])
    #		if(length(indlow) != 0) Pr_exceed[indlow] <- 1
    #		if(length(indhigh) != 0) Pr_exceed[indhigh] <- 0
    # 计算累积分布函数 即分布函�?
    if (distrib == "normal") Pr_exceed <- 1 - ptnorm(int, mean = intensity[i], sd = sigma,
                                                  intensity[i] - nsigma * sigma, intensity[i] + nsigma * sigma)
    rate_mi[i,] <- rate[i] * Pr_exceed
  }
  for (i in 1:nbi)
    rate_all[i] <- sum(rate_mi[, i])
  Pr <- 1 - exp(-rate_all)
  indNA <- which(Pr == 0)
  if (length(indNA) != 0) Pr[indNA] <- NA
  return(data.frame(i = int, rate = rate_all, Pr = Pr))
}

fct_MDG <- function(int, vuln_index) {
  # Mean Damage Grade, Giovinnazzi and Lagomarsino (2004; 2006)
  # Formula SERIANEX, AP5000, pp. 87 with vi-13.1
  # Formula Landtwing, pp. 15 with vi-13.2
  muD <- 5 / 2 * (1 + tanh((int + 6.25 * vuln_index - 13.1) / 2.3))

  # calibration - reduction factor (SERIANEX, AP5000, pp. 97)
  # input: i1, i2
  #	if(int < i1) muD <- muD*0
  #	if(int >= i1 & int <= i2) muD <- muD*interp1(c(i1,i2), c(0,1), int, method="linear")^2
  #	if(int > i2) muD <- muD*1
  return(muD)
}

fct_DGD_binomial <- function(MDG) {
  # Damage Grade Distribution, Braga et al. (1982)
  # Binomial PDF
  dg <- seq(0, 5, 1) #EMS-98 damage grades
  Pr <- numeric(6)
  for (i in 1:6)
    Pr[i] <- factorial(5) / (factorial(dg[i]) * factorial(5 - dg[i])) * ((MDG / 5) ^ dg[i]) * ((1 - MDG / 5) ^ (5 - dg[i]))
  return(Pr)
}

# 计算hazard curve 之后会用到，现在用不�?
V <- 3e4 #30000m3
mmin <- 3
Mmax <- 7
mi <- seq(mmin, Mmax, 0.1) #41个数
nsigma <- 3
a <- c(-2.6, -3.2, -2.0, -1.4, -2.4, -3.8, -3.1, -0.5, -0.9, 0.1, -4.2, -2.8, -1.6)
b <- c(0.7, 0.8, 1.4, 0.9, 1.1, 2.2, 1.8, 1.1, 0.8, 1.6, 1.1, 0.8, 1.0)
nsite <- length(a) # 13
inti <- seq(2, 12, .1)
ni <- length(inti) # 101
nIPE <- 4
#13,4,6,81,101
hazcurve <- array(NA, dim = c(nsite, nIPE, nz, nd, ni))
for (i in 1:nsite)
  #13
  for (j in 1:nIPE)
    #4
    for (k in 1:nz)
      #6
      for (l in 1:nd) {
        #81
        rate <- abs(diff(10 ^ (a[i] - b[i] * mi) * V))
        dhyp <- sqrt(z[k] ^ 2 + d[l] ^ 2)
        if (j == 1) Int <- IPE1(mi, dhyp)
        if (j == 2) Int <- IPE2(mi, dhyp)
        if (j == 3) Int <- IPE3(mi, dhyp)
        if (j == 4) Int <- IPE4(mi, dhyp)
        hazcurve[i, j, k, l,] <- fct_hazardcurve(rate, Int$median, Int$sigma, nsigma, inti, "normal")$Pr
      }

#计算不同建筑等级下的死亡率的不同

Vi <- c(0.4, 0.6, 0.75, 0.9) # class D, C, B, A

nVi <- length(Vi) #4
MDGval <- numeric(nVi * ni);
#4*101
dim(MDGval) <- c(nVi, ni)
DG_pr <- numeric(nVi * ni * 6);
dim(DG_pr) <- c(nVi, ni, 6)
Deathval <- numeric(nVi * ni);
dim(Deathval) <- c(nVi, ni)
# 计算死亡�?
for (i in 1:nVi)
  # nVi <- 4
  for (j in 1:ni) {
    # ni <- 101
    MDGval[i, j] <- fct_MDG(inti[j], Vi[i])
    DG_pr[i, j,] <- fct_DGD_binomial(MDGval[i, j])
    Deathval[i, j] <- 0.00001 * DG_pr[i, j, (2 + 1)]
    + 0.00002 * DG_pr[i, j, (3 + 1)]
    + 0.0002 * DG_pr[i, j, (4 + 1)] + 0.1 * DG_pr[i, j, (5 + 1)]
  }

# 不同的IR下的井失败的概率
safety.Pr <- 1e-6
safety.Pr2 <- 1e-5
deathi <- 10 ^ seq(-9, 0, 0.01)

limit <- 10 ^ (-log10(deathi) + log10(safety.Pr))
limit2 <- 10 ^ (-log10(deathi) + log10(safety.Pr2))

# 井失败的概率
## probability of well failure ##
pr_fail <- array(NA, dim = c(2, nVi, nz, nd))
passfail <- array("green", dim = c(2, nVi, nIPE, nsite, nz, nd))

for (zz in 1:nz) {

  for (xx in 1:nd) {

    for (i in 1:nVi) {

      for (j in 1:nIPE) {

        for (k in 1:nsite) {
          #get p from all sites, max Mmax, all IPEs
          lim.check1 <- 10 ^ (-log10(Deathval[i,]) + log10(safety.Pr))
          indcheck1 <- which(hazcurve[k, j, zz, xx,] >= lim.check1)
          if (length(indcheck1) > 0) passfail[1, i, j, k, zz, xx] <- "red"

          lim.check2 <- 10 ^ (-log10(Deathval[i,]) + log10(safety.Pr2))
          indcheck2 <- which(hazcurve[k, j, zz, xx,] >= lim.check2)
          if (length(indcheck2) > 0) passfail[2, i, j, k, zz, xx] <- "red"
        }
      }
      pr_fail[1, i, zz, xx] <- length(which(passfail[1, i,,, zz, xx] == "red")) / (nsite * nIPE)
      pr_fail[2, i, zz, xx] <- length(which(passfail[2, i,,, zz, xx] == "red")) / (nsite * nIPE)
    }
  }
}

# 更新价格
#按照四种建筑等级进行计算价格
## update price ##
## CPT ##
#mean estimates
lambda <- 2.18
alpha <- 0.78
beta <- 0.82
gamma <- 0.72
delta <- 0.77

# case_norm <- 1 means 1mmt
# case_norm <- 2 means 10mmt
case_norm <- 1
Paverse11 <- array(NA, dim = c(nz, nd)) # nz <- 6; nd <- 81
Paverse21 <- array(NA, dim = c(nz, nd))
Paverse31 <- array(NA, dim = c(nz, nd))
Paverse41 <- array(NA, dim = c(nz, nd))
for (zz in 1:nz) {
  # 6
  for (xx in 1:nd) {
    #81
    E <- E.combL0_z[zz, xx]
    C.TLS <- costs.well_z[zz] + C.frac
    C <- costs.combL0_z[zz, xx]
    pi1 <- pr_fail[case_norm, 1, zz, xx]
    pi2 <- pr_fail[case_norm, 2, zz, xx]
    pi3 <- pr_fail[case_norm, 3, zz, xx]
    pi4 <- pr_fail[case_norm, 4, zz, xx]

    wp1 <- (1 - pi1) ^ gamma / ((1 - pi1) ^ gamma + (1 - (1 - pi1)) ^ gamma) ^ (1 / gamma)
    wp2 <- (1 - pi2) ^ gamma / ((1 - pi2) ^ gamma + (1 - (1 - pi2)) ^ gamma) ^ (1 / gamma)
    wp3 <- (1 - pi3) ^ gamma / ((1 - pi3) ^ gamma + (1 - (1 - pi3)) ^ gamma) ^ (1 / gamma)
    wp4 <- (1 - pi4) ^ gamma / ((1 - pi4) ^ gamma + (1 - (1 - pi4)) ^ gamma) ^ (1 / gamma)
    wm1 <- pi1 ^ delta / (pi1 ^ delta + (1 - pi1) ^ delta) ^ (1 / delta)
    wm2 <- pi2 ^ delta / (pi2 ^ delta + (1 - pi2) ^ delta) ^ (1 / delta)
    wm3 <- pi3 ^ delta / (pi3 ^ delta + (1 - pi3) ^ delta) ^ (1 / delta)
    wm4 <- pi4 ^ delta / (pi4 ^ delta + (1 - pi4) ^ delta) ^ (1 / delta)
    Paverse11[zz, xx] <- 1 / E * (((wm1 * lambda * C.TLS ^ beta) / wp1) ^ (1 / alpha) + C)
    Paverse21[zz, xx] <- 1 / E * (((wm2 * lambda * C.TLS ^ beta) / wp2) ^ (1 / alpha) + C)
    Paverse31[zz, xx] <- 1 / E * (((wm3 * lambda * C.TLS ^ beta) / wp3) ^ (1 / alpha) + C)
    Paverse41[zz, xx] <- 1 / E * (((wm4 * lambda * C.TLS ^ beta) / wp4) ^ (1 / alpha) + C)
  }
}
# npt <- 486
LCOEaverse1_norm1.gg <- rep(NA, npt)
LCOEaverse2_norm1.gg <- rep(NA, npt)
LCOEaverse3_norm1.gg <- rep(NA, npt)
LCOEaverse4_norm1.gg <- rep(NA, npt)
for (pt in 1:npt) {
  j <- grid.ind$i[pt] # 80�? 1�?80
  i <- grid.ind$j[pt] # 80�? 1:6
  LCOEaverse1_norm1.gg[pt] <- Paverse11[i, j]
  LCOEaverse2_norm1.gg[pt] <- Paverse21[i, j]
  LCOEaverse3_norm1.gg[pt] <- Paverse31[i, j]
  LCOEaverse4_norm1.gg[pt] <- Paverse41[i, j]
}
case_norm <- 2
Paverse12 <- array(NA, dim = c(nz, nd))
Paverse22 <- array(NA, dim = c(nz, nd))
Paverse32 <- array(NA, dim = c(nz, nd))
Paverse42 <- array(NA, dim = c(nz, nd))
for (zz in 1:nz) {
  for (xx in 1:nd) {
    E <- E.combL0_z[zz, xx]
    C.TLS <- costs.well_z[zz] + C.frac
    C <- costs.combL0_z[zz, xx]
    pi1 <- pr_fail[case_norm, 1, zz, xx]
    pi2 <- pr_fail[case_norm, 2, zz, xx]
    pi3 <- pr_fail[case_norm, 3, zz, xx]
    pi4 <- pr_fail[case_norm, 4, zz, xx]

    wp1 <- (1 - pi1) ^ gamma / ((1 - pi1) ^ gamma + (1 - (1 - pi1)) ^ gamma) ^ (1 / gamma)
    wp2 <- (1 - pi2) ^ gamma / ((1 - pi2) ^ gamma + (1 - (1 - pi2)) ^ gamma) ^ (1 / gamma)
    wp3 <- (1 - pi3) ^ gamma / ((1 - pi3) ^ gamma + (1 - (1 - pi3)) ^ gamma) ^ (1 / gamma)
    wp4 <- (1 - pi4) ^ gamma / ((1 - pi4) ^ gamma + (1 - (1 - pi4)) ^ gamma) ^ (1 / gamma)
    wm1 <- pi1 ^ delta / (pi1 ^ delta + (1 - pi1) ^ delta) ^ (1 / delta)
    wm2 <- pi2 ^ delta / (pi2 ^ delta + (1 - pi2) ^ delta) ^ (1 / delta)
    wm3 <- pi3 ^ delta / (pi3 ^ delta + (1 - pi3) ^ delta) ^ (1 / delta)
    wm4 <- pi4 ^ delta / (pi4 ^ delta + (1 - pi4) ^ delta) ^ (1 / delta)
    Paverse12[zz, xx] <- 1 / E * (((wm1 * lambda * C.TLS ^ beta) / wp1) ^ (1 / alpha) + C)
    Paverse22[zz, xx] <- 1 / E * (((wm2 * lambda * C.TLS ^ beta) / wp2) ^ (1 / alpha) + C)
    Paverse32[zz, xx] <- 1 / E * (((wm3 * lambda * C.TLS ^ beta) / wp3) ^ (1 / alpha) + C)
    Paverse42[zz, xx] <- 1 / E * (((wm4 * lambda * C.TLS ^ beta) / wp4) ^ (1 / alpha) + C)
  }
}

LCOEaverse1_norm2.gg <- rep(NA, npt)
LCOEaverse2_norm2.gg <- rep(NA, npt)
LCOEaverse3_norm2.gg <- rep(NA, npt)
LCOEaverse4_norm2.gg <- rep(NA, npt)
for (pt in 1:npt) {
  j <- grid.ind$i[pt]
  i <- grid.ind$j[pt]
  LCOEaverse1_norm2.gg[pt] <- Paverse12[i, j]
  LCOEaverse2_norm2.gg[pt] <- Paverse22[i, j]
  LCOEaverse3_norm2.gg[pt] <- Paverse32[i, j]
  LCOEaverse4_norm2.gg[pt] <- Paverse42[i, j]
}
#LCOE 不同IR下的不同建筑等级下的价格计算完成，合成data.frame
LCOE.gg <- data.frame(LCOE.gg,
                      Pa11 = LCOEaverse1_norm1.gg * 1e2,
                      Pa21 = LCOEaverse2_norm1.gg * 1e2,
                      Pa31 = LCOEaverse3_norm1.gg * 1e2,
                      Pa41 = LCOEaverse4_norm1.gg * 1e2,
                      Pa12 = LCOEaverse1_norm2.gg * 1e2,
                      Pa22 = LCOEaverse2_norm2.gg * 1e2,
                      Pa32 = LCOEaverse3_norm2.gg * 1e2,
                      Pa42 = LCOEaverse4_norm2.gg * 1e2)


pmax <- 40
pmin <- 0

z_target <- -6
price_target <- 6


d11min <- min(d[which(Paverse11[z == -z_target,] * 1e2 < price_target)])
d11max <- max(d[which(Paverse11[z == -z_target,] * 1e2 < price_target)])
d21min <- min(d[which(Paverse21[z == -z_target,] * 1e2 < price_target)])
d21max <- max(d[which(Paverse21[z == -z_target,] * 1e2 < price_target)])
d31min <- min(d[which(Paverse31[z == -z_target,] * 1e2 < price_target)])
d31max <- max(d[which(Paverse31[z == -z_target,] * 1e2 < price_target)])
d41min <- min(d[which(Paverse41[z == -z_target,] * 1e2 < price_target)])
d41max <- max(d[which(Paverse41[z == -z_target,] * 1e2 < price_target)])

d12min <- min(d[which(Paverse12[z == -z_target,] * 1e2 < price_target)])
d12max <- max(d[which(Paverse12[z == -z_target,] * 1e2 < price_target)])
d22min <- min(d[which(Paverse22[z == -z_target,] * 1e2 < price_target)])
d22max <- max(d[which(Paverse22[z == -z_target,] * 1e2 < price_target)])
d32min <- min(d[which(Paverse32[z == -z_target,] * 1e2 < price_target)])
d32max <- max(d[which(Paverse32[z == -z_target,] * 1e2 < price_target)])
d42min <- min(d[which(Paverse42[z == -z_target,] * 1e2 < price_target)])
d42max <- max(d[which(Paverse42[z == -z_target,] * 1e2 < price_target)])


pdf(file = paste(wd, "/", figd, "/figNEW_newLCOE_optimalsiting.pdf", sep = ""))
g11 <- ggplot(data = LCOE.gg, aes(x = x, y = -y)) + #先给出坐�?
#创建栅格
geom_raster(aes(fill = Pa11)) +
# 标图不同值得颜色
scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 6, limits = c(pmin, pmax), na.value = "lightgrey") +
#3d曲面�?2d轮廓
geom_contour(aes(z = Pa11), breaks = c(price_target), col = "black") +
# 绘制参考线
geom_hline(yintercept = z_target, lty = "dashed") +
  geom_vline(xintercept = c(d11min, d11max), lty = "dashed") +

  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "LCOE (averse), class D", x = "d (km)", y = "z (km)", fill = expression(P[averse]))

g12 <- ggplot(data = LCOE.gg, aes(x = x, y = -y)) +
  geom_raster(aes(fill = Pa12)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 6, limits = c(pmin, pmax), na.value = "lightgrey") +
  geom_contour(aes(z = Pa12), breaks = c(price_target), col = "black") +
  geom_hline(yintercept = z_target, lty = "dashed") +
  geom_vline(xintercept = c(d12min, d12max), lty = "dashed") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "LCOE (averse), class D", x = "d (km)", y = "z (km)", fill = expression(P[averse]))

g21 <- ggplot(data = LCOE.gg, aes(x = x, y = -y)) +
  geom_raster(aes(fill = Pa21)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 6, limits = c(pmin, pmax), na.value = "lightgrey") +
  geom_contour(aes(z = Pa21), breaks = c(price_target), col = "black") +
  geom_hline(yintercept = z_target, lty = "dashed") +
  geom_vline(xintercept = c(d21min, d21max), lty = "dashed") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "LCOE (averse), class C", x = "d (km)", y = "z (km)", fill = expression(P[averse]))
g22 <- ggplot(data = LCOE.gg, aes(x = x, y = -y)) +
  geom_raster(aes(fill = Pa22)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 6, limits = c(pmin, pmax), na.value = "lightgrey") +
  geom_contour(aes(z = Pa22), breaks = c(price_target), col = "black") +
  geom_hline(yintercept = z_target, lty = "dashed") +
  geom_vline(xintercept = c(d22min, d22max), lty = "dashed") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "LCOE (averse), class C", x = "d (km)", y = "z (km)", fill = expression(P[averse]))

g31 <- ggplot(data = LCOE.gg, aes(x = x, y = -y)) +
  geom_raster(aes(fill = Pa31)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 6, limits = c(pmin, pmax), na.value = "lightgrey") +
  geom_contour(aes(z = Pa31), breaks = c(price_target), col = "black") +
  geom_hline(yintercept = z_target, lty = "dashed") +
  geom_vline(xintercept = c(d31min, d31max), lty = "dashed") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "LCOE (averse), class B", x = "d (km)", y = "z (km)", fill = expression(P[averse]))
g32 <- ggplot(data = LCOE.gg, aes(x = x, y = -y)) +
  geom_raster(aes(fill = Pa32)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 6, limits = c(pmin, pmax), na.value = "lightgrey") +
  geom_contour(aes(z = Pa32), breaks = c(price_target), col = "black") +
  geom_hline(yintercept = z_target, lty = "dashed") +
  geom_vline(xintercept = c(d32min, d32max), lty = "dashed") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "LCOE (averse), class B", x = "d (km)", y = "z (km)", fill = expression(P[averse]))

g41 <- ggplot(data = LCOE.gg, aes(x = x, y = -y)) +
  geom_raster(aes(fill = Pa41)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 6, limits = c(pmin, pmax), na.value = "lightgrey") +
  geom_contour(aes(z = Pa41), breaks = c(price_target), col = "black") +
  geom_hline(yintercept = z_target, lty = "dashed") +
  geom_vline(xintercept = c(d41min, d41max), lty = "dashed") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "LCOE (averse, 1mmt), class A", x = "d (km)", y = "z (km)", fill = expression(P[averse]))
g42 <- ggplot(data = LCOE.gg, aes(x = x, y = -y)) +
  geom_raster(aes(fill = Pa42)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 6, limits = c(pmin, pmax), na.value = "lightgrey") +
  geom_contour(aes(z = Pa42), breaks = c(price_target), col = "black") +
  geom_hline(yintercept = z_target, lty = "dashed") +
  geom_vline(xintercept = c(d42min, d42max), lty = "dashed") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "LCOE (averse, 10mmt), class A", x = "d (km)", y = "z (km)", fill = expression(P[averse]))

grid.arrange(g41, g42, g31, g32, g21, g22, g11, g12, ncol = 2, nrow = 4)
dev.off()

d31best <- d[which(Paverse31[z == -z_target,] * 1e2 == min(Paverse31[z == -z_target,] * 1e2))]
d32best <- d[which(Paverse32[z == -z_target,] * 1e2 == min(Paverse32[z == -z_target,] * 1e2))]




# 综合计算 计算最优位�?

### synthetic exposure ###
#log(rank) <- 10.53-1.005*log(size)    #Gabaix 1999 Zipf law for cities

rank <- seq(10) # 1�? 2�? 3�? 4�? 5�? 6�? 7�? 8�? 9�? 10
A <- 10
# 这是个啥啊？
# size变量意味着什么呢
size <- round(exp(A - log(rank))) # round()四舍五入
# x 0, 200
xmin <- 0
xmax <- 200
# y 0, 100
ymin <- 0
ymax <- 100

# 厂发�? 10MW
plant_power <- 10 #MW
# 用户用电能量
house_power <- 5 #kW

#建筑数量 size之和
nbuilding <- sum(size)
# 村落数量 10
nsettlement <- length(size)
# size * 5*1000/(10*1000000)


nEGS_all <- ceiling(size * house_power * 1e3 / (plant_power * 1e6)) #upper bound

nsim <- 100

target_power <- sum(size * house_power * 1e-3) #MW

reached_power <- array(NA, dim = c(2, nsim))

for (sim in 48:nsim) {
  # �? 48�?100
  x <- numeric(nbuilding) # 64514�?
  y <- numeric(nbuilding) # 64514�?
  # settlementID是i = 1�?10，重复size[i]�?
  settlementID <- unlist(sapply(1:10, function(i) rep(i, size[i])))
  k <- 1
  # 这里的这个过程是怎么算出来的呢？
  for (i in 1:nsettlement) {
    # nsettlement <- 10
    # 均匀分布，随机取
    x[k] <- runif(1, min = xmin, max = xmax) # �? 0,200之间的均匀分布
    y[k] <- runif(1, min = ymin, max = ymax) # �? 0�?100之间的均匀分布
    k <- k + 1
    for (j in 2:size[i]) {
      x[k] <- x[k - 1] + rnorm(1, sd = 0.05)
      y[k] <- y[k - 1] + rnorm(1, sd = 0.05)
      k <- k + 1
    }
  }
  # x, y 是均匀分布加上正态分�?
  xi <- xmin + seq(xmax - xmin);
  nx <- length(xi) # 200
  yi <- ymin + seq(ymax - ymin);
  ny <- length(yi) # 100
  # 这几个变量的维度一�?
  sectorID <- LCOEmap1 <- LCOEmap2 <- dens <- array(NA, dim = c(nx, ny)) # 100�? 200�?

  for (i in 1:nx) {
    # 100
    for (j in 1:ny) {
      # 200
      # indi 返回一系列下标
      indin <- which(x >= xi[i] - .5 & x < xi[i] + .5 & y >= yi[j] - .5 & y < yi[j] + .5) # 返回下标
      dens[i, j] <- length(indin)

      # distance2building 取最�?
      d2building <- sqrt((xi[i] - x) ^ 2 + (yi[j] - y) ^ 2)
      d2building_min <- min(d2building)
      sectorID[i, j] <- settlementID[which(d2building == d2building_min)]

      if (dens[i, j] == 0) {
        LCOEmap1[i, j] <- Paverse31[z == -z_target, d == round(d2building_min)] * 1e2
        LCOEmap2[i, j] <- Paverse32[z == -z_target, d == round(d2building_min)] * 1e2
      }
    }
  }

  x_list <- y_list <- siting_list1 <- siting_list2 <- sector_list <-
    LCOE_list1 <- LCOE_list2 <- nEGS_list <- dens_list <- numeric(nx * ny) # 100 * 200 = 20000
  k <- 1
  for (i in 1:nx) {
    # 100
    for (j in 1:ny) {
      # 200
      x_list[k] <- xi[i]
      y_list[k] <- yi[j]
      sector_list[k] <- sectorID[i, j]
      LCOE_list1[k] <- LCOEmap1[i, j]
      LCOE_list2[k] <- LCOEmap2[i, j]
      dens_list[k] <- dens[i, j]
      # nEGS表示EGS的数�?
      nEGS_list[k] <- ceiling(size[sectorID[i, j]] * house_power * 1e3 / (plant_power * 1e6))
      k <- k + 1
    }
  }
  #把数据写入文�?
  siting.map <- data.frame(x = x_list, y = y_list, dens = dens_list, sector = sector_list, nEGS = nEGS_list,
                           LCOE1 = LCOE_list1, LCOE2 = LCOE_list2)

  write.table(siting.map, paste(wd, "/siting.txt", sep = ""), row.names = F, col.names = T, quote = F)

  # 得到数据之后进行画图
  pmin <- 4
  pmax <- 10
  pdf(paste(wd, "/", figd, "/figNEW_newLCOE_optimalsiting_classB_", sim, ".pdf", sep = ""))
  # A,B是同类型�?
  # 看的是距离和价格的关�?
  gA <- ggplot(data = data.frame(x = d, y = Paverse31[z == -z_target,] * 1e2)) +
    geom_line(aes(x = x, y = y)) +
    geom_vline(xintercept = c(d31min, d31max), lty = "dashed") +
    geom_vline(xintercept = d31best, lty = "dotted") +
    geom_hline(yintercept = price_target, lty = "dashed") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "IR = 1 mmt", x = "d (km)", y = "LCOE (c/kWh)")

  gB <- ggplot(data = data.frame(x = d, y = Paverse32[z == -z_target,] * 1e2)) +
    geom_line(aes(x = x, y = y)) +
    geom_hline(yintercept = price_target, lty = "dashed") +
    geom_vline(xintercept = c(d32min, d32max), lty = "dashed") +
    geom_vline(xintercept = d32best, lty = "dotted") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "IR = 10 mmt", x = "d (km)", y = "LCOE (c/kWh)")
  # C,D是同类型�?
  # 聚类的画�?
  # 画出建筑群的位置
  gC <- ggplot(data = siting.map, aes(x = x, y = y)) +
    geom_raster(aes(fill = LCOE1)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = price_target, limits = c(pmin, pmax), na.value = "black") +
    geom_contour(aes(z = LCOE1), breaks = price_target, colour = "black", lty = "dotted") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "LCOE", x = "x (km)", y = "y (km)")

  gD <- ggplot(data = siting.map, aes(x = x, y = y)) +
    geom_raster(aes(fill = LCOE2)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = price_target, limits = c(pmin, pmax), na.value = "black") +
    geom_contour(aes(z = LCOE2), breaks = price_target, colour = "black", lty = "dotted") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "LCOE", x = "x (km)", y = "y (km)")
  # E,F是同类型�?
  gE <- ggplot(data = siting.map, aes(x = x, y = y)) +
    geom_raster(aes(fill = sector)) +
    scale_color_brewer(palette = "Spectral") +
    theme_minimal() +
  #  theme(legend.position="none") +
  labs(title = "Sectors", x = "x (km)", y = "y (km)")

  gF <- ggplot(data = siting.map, aes(x = x, y = y)) +
    geom_raster(aes(fill = nEGS)) +
    scale_fill_gradient(low = "white", high = "red", na.value = "black", limits = c(0, 12)) +
    theme_minimal() +
  #  theme(legend.position="none") +
  labs(title = "nEGS", x = "x (km)", y = "y (km)")

  grid.arrange(gA, gB, gE, gF, gC, gD, ncol = 2, nrow = 3)
  dev.off()


  # 自由市场方法
  ######################### free market approach ########################
  # 数圆的方�?
  count_circle <- function(X, map, radius) {
    # 数这个圆的数�? #半径也是输入�?
    # nEGS是建立EGS的数�?
    nEGS <- nrow(X)
    # ns是map下sector中的最大�?
    ns <- max(map$sector)

    #how many times in circle per settlement
    count_inriskcircle_persettlement <- numeric(ns)
    for (i in 1:nEGS) {
      #风险圆，得到了圆的数据框
      circle_risk <- circle(c(X$x[i], X$y[i]), radius)
      for (j in 1:ns) {
        #
        indb <- which(map$sector == j & is.na(map$LCOE1) == T)
        settlement_coords <- data.frame(x = map$x[indb], y = map$y[indb])
        indIN <- inpip(settlement_coords, circle_risk)
        if (length(indIN) != 0) count_inriskcircle_persettlement[j] <- count_inriskcircle_persettlement[j] + 1
      }
    }
    indtoorisky <- which(count_inriskcircle_persettlement > 1)
    # 如果这个等于0，就返回“go”，如果不等�?0，就返回“nogo�?
    if (length(indtoorisky) == 0) res <- "go" else res <- "nogo"
    return(res)
  }



  nEGS.max <- sum(nEGS_all) # 38
  indloc <- numeric(nEGS.max) # indloc 38个数
  siting.map1 <- subset(siting.map, LCOE1 < price_target) # 价格1 小于目标价格
  siting.map2 <- subset(siting.map, LCOE2 < price_target) # 价格2 小于目标价格
  radius1 <- d31best # 34
  radius2 <- d32best # 13
  lim1 <- nrow(siting.map1) # 10782
  lim2 <- nrow(siting.map2) # 19634

  #norm 1mmt
  # 在这里运行第一遍的rdm
  lim <- lim1 # 10782
  radius <- radius1 # 34
  indloc_tmp <- which(siting.map1$LCOE1 == min(siting.map1$LCOE1)) #返回下标
  # 向上取整
  rdm <- ceiling(runif(1) * length(indloc_tmp)) # 随机的一个数�?
  indloc <- indloc_tmp[rdm]

  #以上可以得到等与最小值的下标
  # X_EGS写入data_frame x, y
  X_EGS <- data.frame(x = siting.map1$x[indloc], y = siting.map1$y[indloc])
  # EGSsector 是这个下标下的数�?
  EGSsector <- siting.map1$sector[indloc]
  for (i in 2:nEGS.max) {
    # 2:38
    # 进行到哪一步了
    print(paste(i, "/", nEGS.max))
    # 除过这里的值，其余的值小于其余最小的值得下标
    indloc_tmp <- which(siting.map1$LCOE1[-indloc] == min(siting.map1$LCOE1[-indloc]))
    rdm <- ceiling(runif(1) * length(indloc_tmp))
    indloc_tmp <- indloc_tmp[rdm]

    X_EGS_tmp <- rbind(X_EGS, data.frame(x = siting.map1$x[-indloc][indloc_tmp], y = siting.map1$y[-indloc][indloc_tmp]))

    indloc <- c(indloc, seq(nrow(siting.map1))[-indloc][indloc_tmp])
    # 如果输出是nogo 且indloc的长度小于lim，在内部进行计算之后继续判断
    while (count_circle(X_EGS_tmp, siting.map, radius) == "nogo" & length(indloc) < lim) {
      print(length(indloc))
      indloc_tmp <- which(siting.map1$LCOE1[-indloc] == min(siting.map1$LCOE1[-indloc]))
      rdm <- ceiling(runif(1) * length(indloc_tmp))
      indloc_tmp <- indloc_tmp[rdm]
      X_EGS_tmp <- rbind(X_EGS, data.frame(x = siting.map1$x[-indloc][indloc_tmp], y = siting.map1$y[-indloc][indloc_tmp]))
      indloc <- c(indloc, seq(nrow(siting.map1))[-indloc][indloc_tmp])
    }
    # 如果输出为go，则结束，下一次循环即�?
    if (count_circle(X_EGS_tmp, siting.map, radius) == "go") {
      X_EGS <- X_EGS_tmp

      EGSsector <- c(EGSsector, siting.map1$sector[-indloc][indloc_tmp])
      EGSsector_count <- unlist(sapply(1:10, function(i) length(which(EGSsector == i))))
      indmaxout <- which(nEGS_all - EGSsector_count == 0)

      indrem <- c()
      for (j in 1:length(indmaxout))
        indrem <- c(indrem, which(siting.map1$sector == indmaxout[j]))
      indloc <- unique(c(indloc, indrem))
    }
  }
  X_EGS1 <- X_EGS
  EGSsector_count1 <- EGSsector_count
  missingEGS1 <- sum(nEGS_all - EGSsector_count1)
  reached_power[1, sim] <- sum(size * house_power * 1e-3 - (nEGS_all - EGSsector_count1) * 10)




  #norm 10mmt
  lim <- lim2
  radius <- radius2
  indloc_tmp <- which(siting.map2$LCOE2 == min(siting.map2$LCOE2))
  rdm <- ceiling(runif(1) * length(indloc_tmp))
  indloc <- indloc_tmp[rdm]
  X_EGS <- data.frame(x = siting.map2$x[indloc], y = siting.map2$y[indloc])
  EGSsector <- siting.map2$sector[indloc]
  for (i in 2:nEGS.max) {
    print(paste(i, "/", nEGS.max))
    indloc_tmp <- which(siting.map2$LCOE2[-indloc] == min(siting.map2$LCOE2[-indloc]))
    rdm <- ceiling(runif(1) * length(indloc_tmp))
    indloc_tmp <- indloc_tmp[rdm]
    X_EGS_tmp <- rbind(X_EGS, data.frame(x = siting.map2$x[-indloc][indloc_tmp], y = siting.map2$y[-indloc][indloc_tmp]))
    indloc <- c(indloc, seq(nrow(siting.map2))[-indloc][indloc_tmp])
    while (count_circle(X_EGS_tmp, siting.map, radius) == "nogo" & length(indloc) < lim) {
      print(length(indloc))
      indloc_tmp <- which(siting.map2$LCOE2[-indloc] == min(siting.map2$LCOE2[-indloc]))
      rdm <- ceiling(runif(1) * length(indloc_tmp))
      indloc_tmp <- indloc_tmp[rdm]
      X_EGS_tmp <- rbind(X_EGS, data.frame(x = siting.map2$x[-indloc][indloc_tmp], y = siting.map2$y[-indloc][indloc_tmp]))
      indloc <- c(indloc, seq(nrow(siting.map2))[-indloc][indloc_tmp])
    }
    if (count_circle(X_EGS_tmp, siting.map, radius) == "go") {
      X_EGS <- X_EGS_tmp

      EGSsector <- c(EGSsector, siting.map2$sector[-indloc][indloc_tmp])
      EGSsector_count <- unlist(sapply(1:10, function(i) length(which(EGSsector == i))))
      indmaxout <- which(nEGS_all - EGSsector_count == 0)

      indrem <- c()
      for (j in 1:length(indmaxout))
        indrem <- c(indrem, which(siting.map2$sector == indmaxout[j]))
      indloc <- unique(c(indloc, indrem))
    }
  }
  X_EGS2 <- X_EGS
  EGSsector_count2 <- EGSsector_count
  missingEGS2 <- sum(nEGS_all - EGSsector_count2)
  reached_power[2, sim] <- sum(size * house_power * 1e-3 - (nEGS_all - EGSsector_count2) * 10)


  # free market 是为了确定供电的范围
  pdf(paste(wd, "/", figd, "/figNEW_newLCOE_optimalsiting_classB_freemarket_", sim, ".pdf", sep = ""))
  riskinit_zones <- data.frame(x = c(), y = c(), id = c())
  for (i in 1:nrow(X_EGS1)) {
    coords <- circle(c(X_EGS1$x[i], X_EGS1$y[i]), radius1)
    riskinit_zones <- rbind(riskinit_zones, data.frame(x = coords$x, y = coords$y, id = rep(i, 100)))
  }

  gA <- ggplot(data = siting.map, aes(x = x, y = y)) +
    geom_raster(aes(fill = LCOE1)) +
    geom_point(data = X_EGS1, pch = 3) +
    geom_path(data = riskinit_zones, aes(group = id)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = price_target, limits = c(pmin, pmax), na.value = "black") +
    geom_contour(aes(z = LCOE1), breaks = price_target, colour = "black", lty = "dotted") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "IR = 1 mmt", x = "x (km)", y = "y (km)")

  riskinit_zones <- data.frame(x = c(), y = c(), id = c())
  for (i in 1:nrow(X_EGS2)) {
    coords <- circle(c(X_EGS2$x[i], X_EGS2$y[i]), radius2)
    riskinit_zones <- rbind(riskinit_zones, data.frame(x = coords$x, y = coords$y, id = rep(i, 100)))
  }
  gB <- ggplot(data = siting.map, aes(x = x, y = y)) +
    geom_raster(aes(fill = LCOE2)) +
    geom_point(data = X_EGS2, pch = 3) +
    geom_path(data = riskinit_zones, aes(group = id)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = price_target, limits = c(pmin, pmax), na.value = "black") +
    geom_contour(aes(z = LCOE2), breaks = price_target, colour = "black", lty = "dotted") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "IR = 10 mmt", x = "x (km)", y = "y (km)")

  gC <- ggplot(data = data.frame(x = reached_power[1,])) +
    geom_histogram(aes(x), binwidth = 10) +
    geom_vline(xintercept = target_power, col = "darkred", lty = "dashed") +
    xlim(min = 0, max = 350) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Reached power", x = "E (MW)", y = "count")

  gD <- ggplot(data = data.frame(x = reached_power[2,])) +
    geom_histogram(aes(x), binwidth = 10) +
    geom_vline(xintercept = target_power, col = "darkred", lty = "dashed") +
    xlim(min = 0, max = 350) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Reached power", x = "E (MW)", y = "count")

  grid.arrange(gA, gB, gC, gD, ncol = 2, nrow = 3)
  dev.off()


}




