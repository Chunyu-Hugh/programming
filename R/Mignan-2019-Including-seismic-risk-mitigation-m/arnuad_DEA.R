# 差分进化算法

  ######################### differential evolution algo ########################
  ## COST FUNCTION ##
  # 适应度函数
  # 需要知道最后算了啥
  costfunction <- function(X, map, testID, radius_risk) {
    # 维度的第二个
    nEGS <- dim(X)[2]

    #extract settlements
    # 返回下标
    indb <- which(map$sector == testID & is.na(map$LCOE1) == T) 
    # 把这些下标的内容整理成一个数据框
    Xsettlement <- data.frame(x = map$x[indb], y = map$y[indb], dens = map$dens[indb])
    # dens的和
    nb <- sum(Xsettlement$dens)
    # 最大sector
    ns <- max(map$sector)

    #cost parameters
    Power_settlement <- 5e-3 * nb #MW

    #check risk circle building overlap
    # 检查风险圆是否重合
    count_inriskcircle_persettlement <- numeric(ns)
    for (i in 1:nEGS) {
        # 风险圆，返回值是圆的坐标数据
      circle_risk <- circle(X[(1:2), i], radius_risk)

      for (j in 1:ns) { #sector最大值
        indb <- which(map$sector == j & is.na(map$LCOE1) == T)
        settlement_coords <- data.frame(x = map$x[indb], y = map$y[indb])
        indIN <- inpip(settlement_coords, circle_risk)
        if (length(indIN) != 0) count_inriskcircle_persettlement[j] <- count_inriskcircle_persettlement[j] + 1
      }
    }
    indtoorisky <- which(count_inriskcircle_persettlement > 1)
    if (length(indtoorisky) == 0) riskthreshold <- 0 else riskthreshold <- Inf
    
    # 这个公式是怎么来的，而且这个公式下面的值不完整啊
    # 求val的最值
    val <- (Power_settlement - sum(X[3,])) + riskthreshold
    return(val)
    
  }

  ngen <- 10
  safetynorm <- 1 #1=1micromort or 2=10micromort
  testID <- 2 # only for sector 2 for now

  ## PARAMETERS ##
  # 去除唯一性
  nEGS <- unique(siting.map$nEGS[siting.map$sector == testID])


  if (safetynorm == 1) {
    #找到下标
    indloc <- which(siting.map$sector == testID & siting.map$LCOE1 < price_target)
    radius_risk <- d31best
  }
  if (safetynorm == 2) {
    indloc <- which(siting.map$sector == testID & siting.map$LCOE2 < price_target)
    radius_risk <- d32best
  }

  Cr <- 0.9 #crossover probability
  NP <- 10 * nEGS #number of population members
  Fcrit <- sqrt((1 - Cr / 2) / NP) 
  FF <- 0.5
  X_EGS <- array(NA, dim = c(NP, ngen, 3, nEGS))
  cost <- array(NA, dim = c(NP, ngen))

  ### init (gen = 1) ###
  # 第一代
  for (i in 1:NP) {
    rdm <- ceiling(runif(nEGS) * length(indloc))
    X_EGS[i, 1, 1,] <- siting.map$x[indloc][rdm] #x - random uniform in possible space
    X_EGS[i, 1, 2,] <- siting.map$y[indloc][rdm] #y - random uniform in possible space
    X_EGS[i, 1, 3,] <- rep(plant_power, nEGS) #power in MW - best option
    # 进化有问题，维度改变，不知道用的是谁
    # 直接算出适应度函数，看是否符合条件
    cost[i, 1] <- costfunction(X_EGS[1,,,], siting.map, testID, radius_risk)
  }

  ### evolution (gen > 1) ###
  for (g in 2:ngen) {
    # trial populations
    for (i in 1:NP) {
      jrdm <- ceiling(runif(1) * nEGS)
      # 向上取整
      indrdm1 <- ceiling(runif(1) * NP);
      while (indrdm1 == i)
        indrdm1 <- ceiling(runif(1) * NP)
      indrdm2 <- ceiling(runif(1) * NP);
      while (indrdm2 == i | indrdm2 == indrdm1)
        indrdm1 <- ceiling(runif(1) * NP)
      indrdm3 <- ceiling(runif(1) * NP);
      while (indrdm3 == i | indrdm3 == indrdm1 | indrdm3 == indrdm2)
        indrdm1 <- ceiling(runif(1) * NP)

      # generate trial vectors
      for (j in 1:nEGS) {
        # mutation
        # 变异操作
        V <- numeric(3)
        V[1] <- X_EGS[indrdm1, g - 1, 1, j] + FF * (X_EGS[indrdm2, g - 1, 1, j] - X_EGS[indrdm3, g - 1, 1, j]) #x
        V[2] <- X_EGS[indrdm1, g - 1, 2, j] + FF * (X_EGS[indrdm2, g - 1, 2, j] - X_EGS[indrdm3, g - 1, 2, j]) #y
        V[3] <- X_EGS[indrdm1, g - 1, 3, j] + FF * (X_EGS[indrdm2, g - 1, 3, j] - X_EGS[indrdm3, g - 1, 3, j]) #power

        #crossover
        # 交叉
        U <- numeric(3)
        if (runif(1) <= Cr | j == jrdm) U <- V else U <- X_EGS[i, g - 1,, j]
      }
    }
    # select new populations
    for(m in 1:NP){
        cost[m,g] <- costfunction(U[1,,],siting.map, testID, radius_risk)
    }

  }



  ## mutation ##
  Xa_rdm <- array(NA, dim = c(3, nEGS))
  Xa_rdm[1,] <- x_list[indloc][ceiling(runif(nEGS) * length(indloc))]
  Xa_rdm[2,] <- y_list[indloc][ceiling(runif(nEGS) * length(indloc))]
  Xa_rdm[3,] <- c(0, plant_power)[ceiling(runif(nEGS) * 2)]
  Xb_rdm <- array(NA, dim = c(3, nEGS))
  Xb_rdm[1,] <- x_list[indID_safe1][ceiling(runif(nEGS) * length(indID_safe1))]
  Xb_rdm[2,] <- y_list[indID_safe1][ceiling(runif(nEGS) * length(indID_safe1))]
  Xb_rdm[3,] <- c(0, plant_power)[ceiling(runif(nEGS) * 2)]
  V <- array(NA, dim = c(2, maxnEGS))
  V[1,] <- Xinit1_EGS[1,] + FF * (Xb_rdm[1,] - Xc_rdm[1,])
  V[2,] <- Xinit1_EGS[2,] + FF * (Xb_rdm[2,] - Xc_rdm[2,])
  V[3,] <- NULL #cannot follow that equation, whether 0 or plant_power
  # check not outside bound or relocate randomly






  pdf(paste(wd, "/", figd, "/figNEW_newLCOE_optimalsiting_classB_DEalgo.pdf", sep = ""))
  riskinit_zones <- data.frame(x = c(), y = c(), id = c())
  for (i in 1:nEGS) {
    coords <- circle(X_EGS[1,, i], radius_risk)
    riskinit_zones <- rbind(riskinit_zones, data.frame(x = coords$x, y = coords$y, id = rep(i, 100)))
  }

  gA <- ggplot(data = siting.map, aes(x = x, y = y)) +
  geom_raster(aes(fill = LCOE1)) +
  geom_point(data = data.frame(x = X_EGS[1, 1,], y = X_EGS[1, 2,]), pch = 3) +
  geom_path(data = riskinit_zones, aes(group = id)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = price_target, limits = c(pmin, pmax), na.value = "black") +
  geom_contour(aes(z = LCOE1), breaks = price_target, colour = "black", lty = "dotted") +
  geom_contour(aes(z = sector), breaks = testID + 1, colour = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "INIT", x = "x (km)", y = "y (km)")

  gB <- ggplot(data = siting.map, aes(x = x, y = y)) +
  geom_raster(aes(fill = LCOE2)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = price_target, limits = c(pmin, pmax), na.value = "black") +
  geom_contour(aes(z = LCOE2), breaks = price_target, colour = "black", lty = "dotted") +
  #  geom_contour(aes(z=sector), breaks = testID+1, colour = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "INIT", x = "x (km)", y = "y (km)")

  grid.arrange(gA, gB, ncol = 2, nrow = 3)
  dev.off()

  