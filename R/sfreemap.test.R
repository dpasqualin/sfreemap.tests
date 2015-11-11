sfreemap.test.perf <- function(tree_seq
                               , species_seq
                               , q_size_seq
                               , n_sim_seq=c(1)
                               , n_tests=5
                               , parallel=TRUE, serial=TRUE
                               , prog="sfreemap"
                               , message=TRUE, file=NULL) {

    res_size <- length(tree_seq) *
                length(species_seq) *
                length(q_size_seq) *
                length(n_sim_seq) *
                ifelse(isTRUE(parallel) && isTRUE(serial), 2, 1)

    result <- create_result_matrix(res_size)

    r_idx <- 1
    for (t in tree_seq) {
        for (s in species_seq) {
            for (q in q_size_seq) {
                trees <- create_trees(t, s, q)
                for (n in n_sim_seq) {

                    if (isTRUE(serial)) {
                        elapsed <- calc_time(trees, FALSE, prog, n_tests, n)
                        data <- c(t, s, q, elapsed, n, "serial")
                        result[r_idx,] <- data
                        if (isTRUE(message)) {
                            print_info(prog, r_idx, res_size, elapsed, t, s, q, n, "serial")
                        }
                        r_idx <- r_idx + 1
                    }

                    if (isTRUE(parallel)) {
                        elapsed <- calc_time(trees, TRUE, prog, n_tests, n)
                        data <- c(t, s, q, elapsed, n, "parallel")
                        result[r_idx,] <- data
                        if (isTRUE(message)) {
                            print_info(prog, r_idx, res_size, elapsed, t, s, q, n, "parallel")
                        }
                        r_idx <- r_idx + 1
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

print_info <- function(prog, r_idx, res_size, elapsed, t, s, q, n, mode) {
    cat("test", (r_idx), "of", res_size)
    cat(" (prog=", prog
          ,", n_trees=", t
          ,", n_species=", s
          ,", q_size=", q
          ,", n_sim=", n
          ,", mode=", mode
          ,"):", sep="")
    cat(" ", elapsed, "s\n", sep="")
}
