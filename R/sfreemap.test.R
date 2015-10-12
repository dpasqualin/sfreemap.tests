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
