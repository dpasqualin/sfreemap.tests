#!/usr/bin/env Rscript

require(sfreemap.tests)
require(sfreemap)
require(sfreemapc)

# constants
progs = c('sfreemap', 'sfreemapc', 'simmap')
n_tests = 7

run <- function(trees, taxa, q, criteria, nsim=c(1)) {
    for (p in progs) {
        parallel <- ifelse(p == 'simmap', FALSE, TRUE)
        output <- paste(p, criteria, 'txt', sep='.')
        sfreemap.test.perf(trees
                            , taxa
                            , q
                            , nsim
                            , prog=p
                            , parallel=parallel
                            , n_tests=n_tests
                            , file=output)
    }
}

#
# increasing states
#
q <- seq(2,20,2)
run(1, 128, q, 'states')

#
# increasing trees
#
trees <- c(1,seq(2,128,2))
run(trees, 128, 4, 'trees') # serial

#
# increasing taxa
#
taxa <- seq(32,1024,32)
run(1, taxa, 4, 'taxa')

print("DONE")
