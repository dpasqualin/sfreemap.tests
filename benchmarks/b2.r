#!/usr/bin/env Rscript
# performance test for simmap in different number of simulations and states

require(sfreemap.tests)
require(sfreemap)
require(sfreemapc)


n_sim <- seq(20,200,20)
q <- seq(2,10,2)
n_testes <- 5

output <- paste('simmap', 'nsim', 'serial', 'txt', sep='.')
res<-sfreemap.test.perf(1, 128, q, n_sim, prog='simmap', n_tests=5, file=output)
