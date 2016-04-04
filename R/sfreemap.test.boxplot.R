sfreemap.test.boxplot <- function(species=128
                                  , Q_simmap='mcmc'
                                  , Q_sfreemap='mcmc'
                                  , trees=NULL
                                  , n_trees=1
                                  , nsim=25
                                  , n_tests=20
                                  , sample_freq=100
                                  , save_to_file=TRUE
                                  , ... ) {


    pi <- 'equal'
    model <- 'ER'

    if (isTRUE(save_to_file)) {
        dest_dir <- 'tests/'
        if (hasArg(dest_dir)) {
            dest_dir <- list(...)$dest_dir
        }
        dir.create(dest_dir, showWarnings=FALSE)
    }

    # simulate a pure birth tree and the character history
    # use scale=1 to scale the tree to 1
    if (is.null(trees)) {
        tree <- create_trees(n_trees, species, 2)
    }

    hist <- simulation.data(tree)

    metric_values <- c("generation", "diff_lmt", "diff_emr"
                        , "transitions", "time_in_a", "time_in_b")

    if ('phylo' %in% class(tree)) {
        states <- tree$states
    } else {
        states <- tree[[1]]$states
    }

    sfreemap_result <- sfreemap(tree, states, method=Q_sfreemap)

    desc <- summary(sfreemap_result)
    sfreemap_mean <- list(lmt=sum(desc$transitions), emr=desc$dwelling_times[1,sort(unique(states))])
    sfreemap_diff <- sfreemap.diff(hist, sfreemap_mean)

    simmap_result <- matrix(0, nrow=n_tests*nsim
            , ncol=length(metric_values)
            , dimnames=list(1:(n_tests*nsim), metric_values))

    run_parallel <- function(n) {
        mtrees <- make.simmap(tree, states, Q=Q_simmap, pi=pi
                                  , model=model, nsim=nsim
                                  , samplefreq=sample_freq
                                  , message=FALSE)

        return(mtrees)
    }

    all_mtrees <- mclapply(1:n_tests, run_parallel, mc.cores=detectCores())

    for (sim_num in 1:n_tests) {
        mtrees <- all_mtrees[[sim_num]]

        for (i in 1:(length(mtrees)/n_trees)) {
            if ('phylo' %in% class(tree)) {
                v <- simmap.mean(mtrees[[i]])
            } else {
                start <- (i-1) * length(tree) + 1
                end <- start + length(tree) - 1
                v <- simmap.mean(mtrees[start:end])
            }
            diff <- sfreemap.diff(hist, v)
            if (Q_simmap == 'mcmc') {
                step <- i*sample_freq
            } else {
                step <- i
            }
            row <- c(step, diff$lmt, diff$emr, v$lmt, v$emr[1], v$emr[2])
            idx <- ((sim_num-1) * nsim) + i
            simmap_result[as.character(idx),] <- row
        }

    }

    outdir_suffix <- format(Sys.time(), "%Y-%m-%d_%H:%M:%OS")
    outdir_suffix <- paste(outdir_suffix, 'boxplot', sep='_')

    out_dir <- create_out_dir(dest_dir, species, Q_simmap, model, outdir_suffix)
    out_file <- create_out_file(out_dir, 'simmap', nsim)

    create_info_file(out_dir, species=species, Q=Q_simmap, pi=pi, model=model
        , trees=trees, n_trees=n_trees, nsim=nsim, n_tests=n_tests
        , sample_freq=sample_freq)

    plot_boxplot(out_dir, 'boxplot_emr_diff.png', simmap_result, 'diff_emr'
                 , 'Generations', 'Error', sfreemap_diff$emr)

    plot_boxplot(out_dir, 'boxplot_lmt_diff.png', simmap_result, 'diff_lmt'
                 , 'Generations', 'Error', sfreemap_diff$lmt)

    plot_boxplot(out_dir, 'boxplot_lmt.png', simmap_result, 'transitions'
                 , 'Generations', 'Transitions', sfreemap_mean$lmt)

    plot_boxplot(out_dir, 'boxplot_emr_a.png', simmap_result, 'time_in_a'
                 , 'Generations', 'Dwelling time', sfreemap_mean$emr[1])

    plot_boxplot(out_dir, 'boxplot_emr_b.png', simmap_result, 'time_in_b'
                 , 'Generations', 'Dwelling time', sfreemap_mean$emr[2])

    write_to_file(out_file, simmap_result, tree, out_dir, hist)

    return(list(simmap=simmap_result, sfreemap=sfreemap_result, tree=tree))

}
