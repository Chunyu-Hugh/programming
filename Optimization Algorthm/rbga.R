library(genalg)

evaluate <- function(x=c()) {
  returnVal = ((x[1]-7)^2 + (x[2]-77)^2 +(x[3]-777)^2 +(x[4]-7777)^2 +(x[5]-77777)^2)
  returnVal
}


rbga.results = rbga(rep(0.0, times = 5), 
                    rep(99999.9, times = 5), 
                    popSize = 200,
                    iters = 250,
                    evalFunc=evaluate, 
                    verbose=TRUE, 
                    mutationChance=0.05)
summary(rbga.results)

plot(rbga.results)
plot(rbga.results, type="hist")
plot(rbga.results, type="vars")