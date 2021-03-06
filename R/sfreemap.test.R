sfreemap.test.perf <- function(tree_seq
                               , species_seq
                               , q_size_seq
                               , n_sim_seq=c(1)
                               , n_tests=5
                               , omp=c(1)
                               , mc.cores=NULL
                               , parallel=TRUE
                               , serial=TRUE
                               , prog="sfreemap"
                               , message=TRUE
                               , fixed_q=FALSE
                               , estimated_q=TRUE
                               , file=NULL) {

    if (!isTRUE(parallel) && length(mc.cores) > 1) {
        stop('Cannot run in serial with many cores... you seem confused')
    }

    res_size <- length(tree_seq) *
                length(species_seq) *
                length(q_size_seq) *
                length(n_sim_seq) *
                length(omp) *
                length(mc.cores) *
                ifelse(isTRUE(parallel) && isTRUE(serial), 2, 1) *
                ifelse(isTRUE(fixed_q) && isTRUE(estimated_q), 2, 1)

    result <- create_result_matrix(res_size)

    r_idx <- 1
    for (t in tree_seq) {
        for (s in species_seq) {
            for (q in q_size_seq) {
                trees <- create_trees(t, s, q)
                for (o in omp) {
                    for (cores in mc.cores) {
                        for (n in n_sim_seq) {

                            if (any(prog != 'sfreemap' && o > 1,
                                    prog == 'simmap' && isTRUE(parallel))) {
                                result[r_idx,] <- rep(0, ncol(result))
                                if (isTRUE(message)) {
                                    cat ('ignoring run where prog=', prog, ' o=',o,' parallel=', parallel, '\n', sep='')
                                }
                                r_idx <- r_idx + 1
                                next
                            }

                            run <- function(mode, q_type) {
                                q_value <- ifelse(isTRUE(q_type), 'fixed', 'estimated')
                                elapsed <- calc_time(trees, mode, prog,
                                                     n_tests, n, q_type, o,
                                                     cores)
                                mode <- ifelse(isTRUE(mode), 'parallel', 'serial')
                                data <- c(t, s, q, elapsed, n, mode, q_value, o, cores)
                                result[r_idx,] <<- data
                                if (isTRUE(message)) {
                                    print_info(prog, r_idx, res_size,
                                               elapsed, t, s, q, n, mode,
                                               q_value, o, cores)
                                }
                                r_idx <<- r_idx + 1
                            }

                            run_in_mode <- function(mode) {
                                if (isTRUE(fixed_q)) {
                                    run(mode, TRUE)
                                }
                                if (isTRUE(estimated_q)) {
                                    run(mode, FALSE)
                                }
                            }

                            if (isTRUE(serial)) {
                                run_in_mode(FALSE)
                            }

                            if (isTRUE(parallel)) {
                                run_in_mode(TRUE)
                            }
                        }
                    }
                }
            }
        }
    }

    if (!is.null(file)) {
        write_to_file(file, result)
    }

    return(result)
}

print_info <- function(prog, r_idx, res_size, elapsed, t, s, q, n, mode,
                       q_value, omp, cores) {

    cat("test", (r_idx), "of", res_size)
    cat(" (prog=", prog
          ,", n_trees=", t
          ,", n_species=", s
          ,", q_size=", q
          ,", n_sim=", n
          ,", mode=", mode
          ,", q_type=", q_value
          ,", omp=", omp
          ,", cores=", cores
          ,"):", sep="")
    cat(" ", elapsed, "s\n", sep="")
}

calc_time <- function(trees, parallel, prog, n_tests, n_sim, fixed_q, omp,
                      cores, remove_outliers=TRUE) {

    if (inherits(trees, 'phylo')) {
        states <- trees$states
    } else {
        states <- trees[[1]]$states
    }

    doit <- function(expr) {
        t_start <- proc.time()
        expr
        t_end <- proc.time()
        elapsed <- (t_end-t_start)[3]
        return (elapsed)
    }

    # Decide whether to estimate or to use Q from the tree
    if (isTRUE(fixed_q)) {
        if (inherits(trees, 'phylo')) {
            Q <- trees$Q
        } else {
            Q <- trees[[1]]$Q
        }
    } else {
        if (prog == 'sfreemap') {
            Q <- NULL
        } else {
            Q <- 'empirical'
        }
    }

    values <- rep(0, n_tests)

    for (i in 1:n_tests) {
        if (prog == 'sfreemapr') {
            t <- doit(sfreemapr.map(trees, states, Q=Q, parallel=parallel))
        } else if (prog == 'sfreemap') {
            t <- doit(sfreemap(trees, states, Q=Q, method='empirical',
                               type='standard', parallel=parallel, omp=omp,
                               mc.cores=cores))
        } else if (prog == 'simmap') {
            t <- doit(make.simmap(trees, states, Q=Q, nsim=n_sim, message=FALSE))
        } else {
            stop('valid for "prog": (simmap|sfreemap|sfreemapr)')
        }
        values[i] <- t
    }

    if (isTRUE(remove_outliers)) {
        values <- remove_outliers(values)
    }

    t_elapsed <- mean(values)

    return(t_elapsed)
}
