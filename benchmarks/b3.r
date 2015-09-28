#!/usr/bin/env Rscript
# performance test for sfreemapc with openmp

require(sfreemap.tests)
require(sfreemapc)

q <- seq(2,20,2)
n_testes <- 5

output <- paste('sfreemap', 'omp', 'single-tree', 'txt', sep='.')
res<-sfreemap.test.perf(1, 128, q, prog='sfreemapc', n_tests=5, omp=TRUE, file=output)
