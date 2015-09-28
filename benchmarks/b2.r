#!/usr/bin/env Rscript
# performance test for simmap in different number of simulations and states

require(sfreemap.tests)
require(phytools)

n_sim <- c(20,50,100,150,200)
q <- seq(12,20,2)
n_testes <- 5

output <- paste('simmap', 'nsim', 'serial', 'txt', sep='.')
res<-sfreemap.test.perf(1, 128, q, n_sim, prog='simmap', n_tests=5, file=output)
