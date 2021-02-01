library(ggplot2)
f_x <- function(x) {
  return(x ^ 4 + 5 * x ^ 3 + 4 * x ^ 2 - 4 * x + 1)
}
g_x <- function(x) {
  return(2 * (x + 1) ^ 2)
}
x_value <- seq(-4, 1, 0.01)
n_length <- length(x_value)
f_value<- array(NA, dim = n_length)
g_value <- array(NA, dim = n_length)
for (i in 1:n_length) {
  f_value[i] <- f_x(x_value[i])
  g_value[i] <- g_x(x_value[i])

}
data_value_x <- data.frame(x = x_value,f = f_value,g = g_value)
data_f_g <- data.frame(x = f_value,y = g_value)
g1 <- ggplot()+
geom_point(aes(x = -2.96, y = f_x(-2.96)))+
geom_point(aes(x = 0.31, y = 0.3))+
geom_point(aes(x = -1.5, y = 4.19))+
geom_point(aes(x = -1, y = g_x(-1)))+
geom_line(data = data_value_x,aes(x=x, y = f),color = "red")+
geom_line(data = data_value_x,aes(x=x,y = g))+
labs( x = "x", y = "f(x) and g(x)")

g2 <- ggplot()+
geom_point(data = data_f_g,aes(x=x, y = y),color = "red")+
labs( x = "f(x)", y = "g(x)")