#!/usr/bin/env Rscript

require(sfreemap.tests)
require(sfreemap)
require(sfreemapc)

# constants
progs = c('sfreemap', 'sfreemapc', 'simmap')
n_tests = 5

run <- function(trees, taxa, q, parallel, criteria) {
    type <- ifelse(isTRUE(parallel), 'parallel', 'serial')
    for (p in progs) {
        # simmap doesn't have a parallel option
        if (p == 'simmap' & isTRUE(parallel)) {
            next
        }
        output <- paste(p, criteria, type, 'txt', sep='.')
        sfreemap.test.perf(trees, taxa, q,
                            , parallel=parallel
                            , n_tests=n_tests
                            , file=output)
    }
}

#
# increasing states
#
trees <- 1
q <- seq(2,20,2)
taxa <- 128
run(trees, taxa, q, FALSE, 'states')

#
# increasing trees
#
trees <- c(1,seq(2,90,2))
q <- 4
taxa <- 128
run(trees, taxa, q, FALSE, 'trees') # serial
run(trees, taxa, q, TRUE, 'trees')  # parallel

#
# increasing taxa
#
trees <- 1
q <- 4
taxa <- seq(32,1024,32)
run(trees, taxa, q, FALSE, 'taxa')
