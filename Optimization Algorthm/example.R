#»ùÒòËã·¨



f<-function(x){
  return (x^2 +5*x+4)
}

library(genalg)

m <- mcga(popsize = 200
          , chsize = 1
          ,crossprob = 1.0
          , mutateprob = 0.01
          ,minval = -5
          , maxval = 5
          , maxiter = 10
          , evalFunc = f)
cat("Best chromosome:\n")
print(m$population[1,])
cat("Cost: ",m$costs[1],"\n")