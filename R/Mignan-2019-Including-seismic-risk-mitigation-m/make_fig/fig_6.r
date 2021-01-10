library(ggplot2)
x <- c(-2.6, -3.2, -2.0, -1.4, -2.4, -3.8, - 3.1, -0.5, -0.9, 0.1, -4.2, -2.8, -1.6)
y <- c(0.7,  0.8,  1.4, 0.9, 1.1, 2.2, 1.8, 1.1, 0.8, 1.6, 1.1, 0.8 ,1.0)
ggplot(data = data.frame(a = x, b = y), aes(x=a, y=b))+
geom_point()+
labs(title = "Underground feedback", x = expression(a(m^-3)))+
theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
