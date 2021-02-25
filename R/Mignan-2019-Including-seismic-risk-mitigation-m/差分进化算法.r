# 这是完整的算法，没有画图而已
cost_function <- function(x) {
  return(3 * cos(x[1] * x[2]) + x[1] + x[2])
}

NP <- 20 #种群数量
D <- 2 #变量的维度
G <- 100 #最大进化代数
f <- 0.5 # 变异因子
CR <- 0.1 # 交叉因子

Xs <- 4 # 上限
Xx <- -4 #下限
# 初始化种群
# 赋初值
x <- array(0, dim = c(D, NP)) # 初始种群
V <- array(0, dim = c(D, NP)) # 变异种群
U <- array(0, dim = c(D, NP)) # 选择种群
final <- array(0, dim = G)
# 赋予初值
x[,] = runif(D * NP, Xx, Xs)
obj_1 <- array(0, dim = NP)
obj <- array(0, dim = NP)
for (m in 1:NP) {
  para <- c()
  para <- append(para, x[, m])
  obj_1[m] <- cost_function(para)
}


final[1] <- min(obj_1)


for (gen in 1:G) {
  for (i in 1:NP) {

    # 变异操作
    r1 <- ceiling(runif(1,0,1) * NP);
    while (r1 == i)
      r1 <- ceiling(runif(1,0,1) * NP)

    r2 <- ceiling(runif(1,0,1) * NP);
    while (r2 == i | r2 == r1)
      r2 <- ceiling(runif(1,0,1) * NP)

    r3 <- ceiling(runif(1,0,1) * NP);
    while (r3 == i | r3 == r1 | r3 == r2)
      r3 <- ceiling(runif(1,0,1) * NP)

    V[, i] <- x[, r1] + f * (x[, r2] - x[, r3])
  }
  #交叉操作
  r <- ceiling(runif(1,0,1) * D)
  for (n in 1:D) {
    #D是变量的维度
    cr <- runif(1)
    if (cr < CR | n == r) U[n,] = V[n,]
    else U[n,] = x[n,]
  }
  #边界处理
  for (n in 1:D) {
    for (m in 1:NP) {
      if (U[n, m] < Xx) {
        U[n, m] <- Xx
      }
      if (U[n, m] > Xs) {
        U[n, m] <- Xs
      }
    }
  }
  # 选择操作
  for (m in 1:NP) {
    para <- c()
    para <- append(para, U[, m])
    obj[m] <- cost_function(para)
  }
  for (m in 1:NP) {
    if (obj[m] < obj_1[m]) {
      x[, m] <- U[, m]
    }
  }
  #计算选择之后的值
  for (m in 1:NP) {
    para <- c()
    para <- append(para, x[, m])
    obj_1[m] <- cost_function(para)
  }
  final[gen + 1] <- min(obj_1)
}