#!/usr/bin/env Rscript

require(sfreemap.tests)
require(sfreemap)
require(sfreemapc)

# constants
n_tests <- 5
nsim <- c(1,10,100)

run <- function(trees, taxa, q, criteria, nsim=c(1)
                     , parallel = TRUE
                     , progs = c('sfreemap', 'sfreemapc', 'simmap')) {
    for (p in progs) {
        parallel <- ifelse(p == 'simmap', FALSE, parallel)
        nsim <- ifelse(p == 'simmap', nsim, c(1))
        output <- paste(p, criteria, 'txt', sep='.')
        sfreemap.test.perf(trees
                            , taxa
                            , q
                            , nsim
                            , prog=p
                            , parallel=parallel
                            , fixed_q=TRUE
                            , n_tests=n_tests
                            , file=output)
    }
}

#
# increasing states
#
#q <- seq(2,20,2)
#run(1, 256, q, 'states', parallel=FALSE, nsim=nsim)

#
# increasing trees
#
trees <- c(1,seq(2,128,2))
run(trees, 256, 4, 'trees', nsim=nsim)

#
# increasing taxa
#
taxa <- seq(32,1024,32)
run(1, taxa, 4, 'taxa', nsim=nsim)

print("DONE")
