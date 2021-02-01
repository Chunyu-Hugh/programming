
# Day24-TSP
library(GA)

# fitness vector
A<-c(0,25,25,10,35)
B<-c(25,0,25,35,25)
C<-c(25,25,0,30,15)
D<-c(10,35,30,0,20)
E<-c(35,25,15,20,0)
city<-rbind(A,B,C,D,E)
colnames(city)<-c('A','B','C','D','E')
data<-as.matrix(city)

#fitness function
tspfitness  <- function(tour, ...) 1/tourLength(tour, ...)

tourLength <- function(tour, distMatrix) {
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[, 2:1]
  sum(route)
}

#GA
GA.fit <- ga(type = "permutation", 
             fitness = tspfitness, 
             distMatrix = data,
             lower  = 1, 
             upper = 5,
             popSize = 100,
             maxiter = 200,
             run = 10, 
             pmutation = 0.8, 
             monitor = NULL)

# Result
summary(GA.fit)

