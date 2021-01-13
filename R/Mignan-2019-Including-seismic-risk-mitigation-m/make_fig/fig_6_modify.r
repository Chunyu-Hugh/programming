list(rm = ls())
library(ggplot2)

x1 <- c(-2.6, -3.2, -2.0, -3.8, -3.1, -0.5, -0.9, -4.2)
y1 <- c(0.7, 0.8, 1.4, 2.2, 1.8, 1.1, 0.8, 1.1)
t1 <- c(rep("D&S13", length(x1)))
d1 <- data.frame(a = x1, b = y1, t = t1)

x2 <- c(-1.4, -2.4, 0.1, -2.8, -1.6)
y2 <- c(0.9, 1.1, 1.6, 0.8, 1.0)
t2 <- c(rep("M&al17", length(x2)))
d2 <- data.frame(a = x2, b = y2, t = t2)

datax <- merge(d1, d2, all = TRUE)
g1 <- ggplot(data = datax, aes(x = a, y = b, group = t, col = t)) +
  geom_point(size = 2) +
  labs(
    title = "(a) Underground feedback",
    x = expression(a(m^-3)), color = "Reference"
  ) +
  theme(plot.title = element_text(hjust = 0.5))



# California
IPE1 <- function(m, r) {
  # Atkinson and Wald, 2007 - California
  c1 <- 12.27
  c2 <- 2.270
  c3 <- 0.1304
  c4 <- -1.30
  c5 <- -0.0007070
  c6 <- 1.95
  c7 <- -0.577
  Rt <- 30
  sigma <- 0.4
  B <- log10(r / Rt)
  int <- c1 + c2 * (m - 6) + c3 * (m - 6)^2 + c4 * log10(r) + c5 * r + c6 * B + c7 * m * log10(r)
  return(c(int, sigma))
}
# Central/Eastern U.S.A
IPE2 <- function(m, r) {
  # Atkinson and Wald, 2007 - Central/Eastern US
  c1 <- 11.72
  c2 <- 2.36
  c3 <- 0.1155
  c4 <- -0.44
  c5 <- -0.002044
  c6 <- 2.31
  c7 <- -0.479
  Rt <- 80
  sigma <- 0.4
  B <- log10(r / Rt)
  int <- c1 + c2 * (m - 6) + c3 * (m - 6)^2 + c4 * log10(r) + c5 * r + c6 * B + c7 * m * log10(r)
  return(c(int, sigma))
}
# Switzerland
IPE3 <- function(m, d) {
  # m moment magnitude, d hypocentral distance (km)
  # int observed Mercalli Intensity
  a <- -0.69182
  b <- 0.00084
  alpha <- 0.7317
  beta <- 1.2567
  c0 <- beta
  c1 <- alpha
  c2 <- -alpha * a
  c3 <- -alpha * b
  int <- (m - c2 * log(d / 30) - c3 * (d - 30) - c0) / c1
  sigma <- 1
  return(c(int, sigma))
}
# Global
IPE4 <- function(m, d) {
  # Allen et al. (2012), pp. 418, Eqs. 2 and 3 (here without site effect), pp. 423, Eq. 6
  # m moment magnitude, d hypocentral distance (km)
  # int generic Cancani scale (12 degrees)
  c0 <- 2.085
  c1 <- 1.428
  c2 <- -1.402
  c4 <- 0.078
  m1 <- -0.209
  m2 <- 2.042
  Rm <- m1 + m2 * exp(m - 5)

  # 	if(d <= 50) int <- c0+c1*m+c2*log(sqrt(d^2+Rm^2))
  # 	else int <- c0+c1*m+c2*log(sqrt(d^2+Rm^2))+c4*log(d/50)

  int <- c0 + c1 * m + c2 * log(sqrt(d^2 + Rm^2)) + c4 * log(d / 50)

  s1 <- 0.82
  s2 <- 0.37
  s3 <- 22.9
  sigma <- s1 + s2 / (1 + (d / s3)^2)
  return(c(int, sigma))
}

# 需要循环进行取
# 输入参数
dis <- seq(0, 80, 0.01)
nDis <- length(dis)
I <- array(NA, dim = c(4, nDis))
re <- c(0, 0)
si <- c(rep(0, 4))
for (i in 1:4) {
  for (j in 1:nDis) {
    if (i == 1) {
      re <- IPE1(6.5, sqrt(4^2 + dis[j]^2))
      I[i, j] <- re[1]
      if (j == nDis) {
        si[i] <- re[2]
      }
    }
    if (i == 2) {
      re <- IPE2(7, sqrt(4^2 + dis[j]^2))
      I[i, j] <- re[1]
      if (j == nDis) {
        si[i] <- re[2]
      }
    }
    if (i == 3) {
      re <- IPE3(6, sqrt(4^2 + dis[j]^2))
      I[i, j] <- re[1]
      if (j == nDis) {
        si[i] <- re[2]
      }
    }
    if (i == 4) {
      re <- IPE4(6, sqrt(4^2 + dis[j]^2))
      I[i, j] <- re[1]
      if (j == nDis) {
        si[i] <- re[2]
      }
    }
  }
}

type <- c(rep("California", nDis), rep("Central/Eastern U.S.A", nDis), rep("Switzerland", nDis), rep("Global", nDis))
cI <- c(I[1, ], I[2, ], I[3, ], I[4, ])
data_Shaking <- data.frame(d = seq(0, 80, 0.01), cI = cI, type = type)
g2 <- ggplot() +
  geom_line(data = data_Shaking, aes(x = d, y = cI, group = type, col = type)) +
  #geom_ribbon(aes(x = seq(0, 80, 0.01),ymin = I[1, ] * (1 - si[1]), ymax = I[1, ] * (1 + si[1]),type = "California"),  alpha = .3, fill = "orange")




g2 <- ggplot(data = data.frame(x = seq(0, 80, 0.01), y = I[1, ]), aes(x = x, y = y)) +
  geom_line(col = "green") +
  geom_ribbon(aes(ymin = I[1, ] * (1 - si[1]), ymax = I[1, ] * (1 + si[1])), alpha = .3, fill = "orange") +

  # geom_line(data = data.frame(x = seq(0, 80, 0.01), y = I[2, ]), aes(x = x, y = y), col = "darkred") +
  geom_ribbon(aes(ymin = I[2, ] * (1 - si[2]), ymax = I[2, ] * (1 + si[2])), alpha = .3, fill = "#00f7ff") +

  # geom_line(data = data.frame(x = seq(0, 80, 0.01), y = I[3, ]), aes(x = x, y = y), col = "blue") +
  geom_ribbon(aes(ymin = I[3, ] * (1 - si[3]), ymax = I[3, ] * (1 + si[3])), alpha = .3, fill = "#09e44b") +

  # geom_line(data = data.frame(x = seq(0, 80, 0.01), y = I[4, ]), aes(x = x, y = y), col = "black") +
  geom_ribbon(aes(ymin = I[4, ] * (1 - si[4]), ymax = I[4, ] * (1 + si[4])), alpha = .4, fill = "#C8FF00") +
  labs(title = "(c) Shaking intensity", x = "d (km)", y = "I", col = "Region") +
  ylim(min = 0, max = 10) +
  theme(plot.title = element_text(hjust = 0.5))

mean_I <- c(rep(0, nDis))
for (i in 1:nDis) {
  mean_I[i] <- mean(I[, i])
}
muD1 <- (1 + tanh((mean_I + 6.25 * .9 - 13.1) / 2.3)) * 2.5
muD2 <- (1 + tanh((mean_I + 6.25 * .75 - 13.1) / 2.3)) * 2.5
muD3 <- (1 + tanh((mean_I + 6.25 * .6 - 13.1) / 2.3)) * 2.5
muD4 <- (1 + tanh((mean_I + 6.25 * .4 - 13.1) / 2.3)) * 2.5
data <- data.frame(I = rep(mean_I, 4), D = c(muD1, muD2, muD3, muD4), Building_class = c(rep("A", nDis), rep("B", nDis), rep("C", nDis), rep("D", nDis)))
g3 <- ggplot(data) +
  geom_smooth(aes(x = I, y = D, group = Building_class, col = Building_class))
# geom_line(data = data.frame(x = mean_I, y = muD2), aes(x = x,y =y),col = "red")+
# geom_line(data = data.frame(x = mean_I, y = muD3), aes(x = x,y =y),col = "green")+
# geom_line(data = data.frame(x = mean_I, y = muD4), aes(x = x,y =y),col = "yellow")+
# scale_colour_discrete(breaks = c('A','B','C','D'), labels = c('W','X','Y','Z'))