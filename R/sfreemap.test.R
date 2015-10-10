sfreemap.test.perf <- function(tree_seq, species_seq, q_size_seq
                               , n_sim_seq=c(1), omp=FALSE
                               , n_tests=5, parallel=TRUE, prog='sfreemap'
                               , message=TRUE, file=NULL) {

    res_size <- length(tree_seq) * length(species_seq) * length(q_size_seq) * length(n_sim_seq)

    result <- create_result_matrix(res_size)

    r_idx <- 0
    for (t in tree_seq) {
        for (s in species_seq) {
            for (q in q_size_seq) {
                trees <- create_trees(t, s, q)
                for (n in n_sim_seq) {
                    if (isTRUE(message)) {
                        cat('test', (r_idx+1), 'of', res_size)
                        cat(' (n_trees=', t
                              ,', n_species=', s
                              ,', q_size=', q
                              ,', n_sim=', n
                              ,'): ', sep='')
                    }
                    elapsed <- calc_time(trees, parallel, prog, n_tests, n, omp)
                    data <- c(t, s, q, elapsed, n)
                    r_idx <- r_idx + 1
                    result[r_idx,] <- data
                    if (isTRUE(message)) {
                        cat (elapsed, 's\n', sep='')
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

simmap.plot_tree <- function(trees, nsim, samplefreq=100) {
    if (is.null(trees[[1]]$maps)) {
        # calculate first and than plot
    }

    total <- length(trees)/nsim
    mean_trees <- list()
    for (i in 1:total) {
        t <- trees[[i]]
        range <- ((i-1)*nsim+1) : (i*nsim)
        # remember that t$maps will be wrong here, but it's ok because we
        # don't need it here
        mapped.edge <- lapply(trees[range], function(x) x$mapped.edge)
        t$mapped.edge <- Reduce('+', mapped.edge) / length(mapped.edge)
        mean_trees[[i]] <- t
    }
    class(mean_trees) <- 'multiPhylo'

    base_tree <- mean_trees[[1]]
    return(mean_trees)
    for (state in colnames(base_tree$mapped.edge)) {
        cat('plotting for state ', state, '\n')
        x11()
        sfreemapc::sfreemap.plot_tree(base_tree, mean_trees, state)
    }
}
