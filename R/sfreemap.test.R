sfreemap.test.boxplot <- function(species=100
										   , Q='mcmc'
										   , pi='equal'
										   , model='ER'
										   , trees=NULL
                                           , n_trees=1
										   , nsim=50
										   , n_tests=5
										   , sample_freq=100
										   , save_to_file=TRUE
						  				   , ... ) {


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

	simmap_result <- matrix(0, nrow=n_tests*nsim
			, ncol=length(metric_values)
			, dimnames=list(1:(n_tests*nsim), metric_values))

    if (class(tree) == 'phylo') {
        states <- tree$states
    } else {
        states <- tree[[1]]$states
    }

	for (sim_num in 1:n_tests) {
		mtrees <- make.simmap(tree, states, Q=Q, pi=pi
								  , model=model, nsim=nsim
								  , samplefreq=sample_freq
								  , message=FALSE)

		for (i in 1:(length(mtrees)/n_trees)) {
            if (class(tree) == 'phylo') {
			    v <- simmap.mean(mtrees[[i]])
            } else {
                start <- (i-1) * length(tree) + 1
                end <- start + length(tree) - 1
                v <- simmap.mean(mtrees[start:end])
            }
			diff <- sfreemap.diff(hist, v)
			row <- c(i*sample_freq, diff$lmt, diff$emr, v$lmt, v$emr[1], v$emr[2])
			idx <- ((sim_num-1) * nsim) + i
			simmap_result[as.character(idx),] <- row
		}

	}

	sfreemap_result <- sfreemap.map(tree, states, Q=Q
									, n_simulations=nsim
									, sample_freq=sample_freq)

	desc <- sfreemap.describe(sfreemap_result)
	sfreemap_mean <- list(lmt=sum(desc$transitions), emr=desc$dwelling_times)
	sfreemap_diff <- sfreemap.diff(hist, sfreemap_mean)

	outdir_suffix <- format(Sys.time(), "%Y-%m-%d_%H:%M:%OS")
	outdir_suffix <- paste(outdir_suffix, 'boxplot', sep='_')

	out_dir <- create_out_dir(dest_dir, species, 'mcmc', 'ER', outdir_suffix)
	out_file <- create_out_file(out_dir, 'simmap', nsim)

	create_info_file(out_dir, species=species, Q=Q, pi=pi, model=model
		, trees=trees, n_trees=n_trees, nsim=nsim, n_tests=n_tests
		, sample_freq=sample_freq)

	plot_boxplot(out_dir, 'boxplot_emr_diff.png', diff_emr ~ generation
			     , data=simmap_result, 'Generations'
				 , 'Error', 'Error expected markov rewards'
				 , sfreemap_diff$emr)

	plot_boxplot(out_dir, 'boxplot_lmt_diff.png', diff_lmt ~ generation
	     	     , data=simmap_result, 'Generations'
				 , 'Error', 'Error labelled markov transitions'
				 , sfreemap_diff$lmt)

	plot_boxplot(out_dir, 'boxplot_lmt.png', transitions ~ generation
	     	     , data=simmap_result, 'Generations'
				 , 'Transitions', 'Labelled markov transitions'
				 , sfreemap_mean$lmt)

	plot_boxplot(out_dir, 'boxplot_emr_a.png', time_in_a ~ generation
	     	     , data=simmap_result, 'Generations'
				 , 'Dwelling time', 'Expected dwelling time in state "a"'
				 , sfreemap_mean$emr[1])

	plot_boxplot(out_dir, 'boxplot_emr_b.png', time_in_b ~ generation
	     	     , data=simmap_result, 'Generations'
				 , 'Dwelling time', 'Expected dwelling time in state "b"'
				 , sfreemap_mean$emr[2])

	write_to_file(out_file, simmap_result, tree, out_dir, hist)

	return(list(simmap=simmap_result, sfreemap=sfreemap_result, tree=tree))

}

sfreemap.test.perf <- function(tree_seq, species_seq, q_size_seq
								, n_tests=5, parallel=TRUE, prog='sfreemap'
								, message=TRUE, file=NULL) {

	res_size <- length(tree_seq) * length(species_seq) * length(q_size_seq)

	result <- create_result_matrix(res_size)

	r_idx <- 0
	for (t in tree_seq) {
		for (s in species_seq) {
			for (q in q_size_seq) {
				if (isTRUE(message)) {
					cat('test', (r_idx+1), 'of', res_size)
					cat(' (n_trees=', t
						  ,', n_species=', s
						  ,', q_size=', q
						  ,'): ', sep='')
				}
				trees <- create_trees(t, s, q)
				elapsed <- calc_time(trees, parallel, prog, n_tests)
				data <- c(t, s, q, elapsed)
				r_idx <- r_idx + 1
				result[r_idx,] <- data
				if (isTRUE(message)) {
					cat (elapsed, 's\n', sep='')
				}
			}
		}
	}

	if (!is.null(file)) {
		write_to_file(file, result)
	}

	return(result)
}
