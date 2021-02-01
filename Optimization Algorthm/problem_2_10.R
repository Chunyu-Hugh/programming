library(ggplot2)
f_1 <- function(x_1, x_2) {
  return(x_1 ^ 2 + x_2)
}

f_2 <- function(x_1, x_2) {
  return(x_2 ^ 2 + x_1)
}
x_1 <- seq(-10, 10, 1)
x_2 <- seq(-10, 10, 1)
nlength <- length(x_1)

f_1_value <- array(NA, dim = c(nlength, nlength))
f_2_value <- array(NA, dim = c(nlength, nlength))

for (i in 1:nlength) {
  for (j in 1:nlength) {
    f_1_value[i, j] <- f_1(x_1[i], x_2[j])
    f_2_value[i, j] <- f_2(x_1[i], x_2[j])

  }
}

f1 <- c()
f2 <- c()
for (i in 1:nlength) {
  f1 <- c(f1, f_1_value[i,])
  f2 <- c(f2, f_2_value[i,])
}
data <- data.frame(f1 = f1, f2 = f2)
ggplot()+
geom_line(data = data,aes(x=f1, y = f2),color = "red")+
labs( x = "f1(x)", y = "f2(x)")