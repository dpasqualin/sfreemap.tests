sfreemap.test.box_and_whiskers <- function(species=100
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
		# Q used in the simulation
		QS <- matrix(c(-1,1,1,-1), 2, 2)
		rownames(QS)<-colnames(QS)<-letters[1:nrow(QS)]

		topologies <- pbtree(n=species, nsim=n_trees, scale=1)
        if (n_trees > 1) {
			tree <- lapply(topologies, sim.history, Q=QS, message=FALSE)
            class(tree) <- 'multiPhylo'
        } else {
            tree <- sim.history(tree, QS, message=FALSE)
	    }
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
            print(idx)
            print(row)
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

	write_to_file(tree, out_dir, out_file, hist, simmap_result)

	return(list(simmap=simmap_result, sfreemap=sfreemap_result, tree=tree))

}

# print all arguments to a file
create_info_file <- function(out_dir, ...) {
	args <- list(...)
	out_file <- paste(out_dir, 'arguments.txt', sep='/')
	for (i in names(args)) {
		if (nchar(i) > 0) {
			txt <- paste(i, '=', args[[i]])
			write(txt, file=out_file, append=TRUE)
		}
	}
}

plot_boxplot <- function(out_dir, out_file, formula, data, xlab, ylab, title, line_data) {

	output <- paste(out_dir, out_file, sep='/')
	png(output, width=1280, height=720)
	boxplot(formula, data=data, xlab=xlab, ylab=ylab, las=2, main=title)
	abline(h=line_data, col='red')
	dev.off()

}

sfreemap.test <- function(species=100, nsim_values=NULL, Q=NULL, QS=NULL
			     		  , model='ER', pi="estimated", message=FALSE
			     		  , test_sfreemap=TRUE, test_simmap=TRUE
				 		  , save_to_file=TRUE
						  , ... ) {

    if (isTRUE(save_to_file)) {
		dest_dir <- 'tests/'
	    if (hasArg(dest_dir)) {
	        dest_dir <- list(...)$dest_dir
	    }
    	dir.create(dest_dir, showWarnings=FALSE)
    }

	# when Q=mcmc, how often should we sample from the simulation?
	sample_freq <- 100
	if (hasArg(sample_freq)) {
		sample_freq <- list(...)$sample_freq
	}

	n_tests <- 1
	if (hasArg(n_tests)) {
		n_tests <- list(...)$n_tests
	}

	if (is.null(nsim_values)) {
		# powers of two up to phylo
		nsim_values <- Filter(function(x) log2(x)%%1==0, seq(2,2048))
	}

    if (is.null(Q)) {
        # set a default Q matrix
    	Q <- "empirical"
    }

    if (is.null(QS)) {
        # set a default Q matrix
    	QS <- matrix(c(-1,1,1,-1), 2, 2)
    }

    if (is.null(model)) {
    	# 'ER' stands for equal rate. This model is related to the
    	# definition of the Q matrix
        model <- 'ER'
    }

    if (is.matrix(Q) && is.null(rownames(Q))) {
        # set rownames and colnames for Q
	    rownames(Q) <- colnames(Q) <- letters[1:nrow(Q)]
    }

    # simulate a pure birth tree and the character history
    # use scale=1 to scale the tree to 1
	tree <- pbtree(n=species, scale=1)
	tree <- sim.history(tree, QS, message=FALSE)

 	hist <- simulation.data(tree)

	# init result
	res_sfreemap <- NULL
	res_simmap <- NULL

	outdir_suffix <- format(Sys.time(), "%Y-%m-%d_%H:%M:%OS")
	out_dir <- create_out_dir(dest_dir, species, Q, model, outdir_suffix)

	for (sim_num in 1:n_tests) {
		if (isTRUE(test_sfreemap)) {
			if (isTRUE(message)) {
				cat('Running sfreemap..\n')
			}
			res_sfreemap <- sfreemap.test.sfreemap(tree, Q, nsim_values, hist)
			if (isTRUE(save_to_file)) {
				out_file <- create_out_file(out_dir, 'sfreemap', sim_num)
				write_to_file(tree, out_dir, out_file, hist, res_sfreemap)
			}
		}

		if (isTRUE(test_simmap)) {
			if (isTRUE(message)) {
				cat('Running simmap..\n')
			}
			res_simmap <- sfreemap.test.simmap(tree, Q, pi, model, nsim_values, hist, sample_freq)
			if (isTRUE(save_to_file)) {
				out_file <- create_out_file(out_dir, 'simmap', sim_num)
				write_to_file(tree, out_dir, out_file, hist, res_simmap)
			}
		}
	}

	return (list(sfreemap=res_sfreemap, simmap=res_simmap))
}

create_out_dir <- function(dest_dir, species, Q, model, outdir_suffix) {
	q_txt <- ifelse(is.matrix(Q), 'matrix', Q)
	out_dir <- paste(species, q_txt, model, outdir_suffix, sep='_')
	out_dir <- paste(dest_dir, out_dir, sep='/')
	dir.create(out_dir, showWarnings=FALSE)
	return (out_dir)
}

create_out_file <- function(out_dir, method, sim_num) {
	out_file <- paste(method, '_', sim_num, '.txt', sep='')
	out_file <- paste(out_dir, out_file, sep='/')
	return (out_file)
}

save_tree_file <- function(out_dir, tree) {
	out_tree_file <- paste(out_dir, '/tree.nexus', sep='')
	write.tree(tree, file=out_tree_file)
	return (out_tree_file)
}

write_to_file <- function(tree, out_dir, out_file, hist, result) {

	save_tree_file(out_dir, tree)

	txt <- paste('# sim.history ', 0, hist$lmt, paste(hist$emr, collapse=' '))
	write(txt, file=out_file)

	txt <- paste('#',paste(colnames(result), collapse=','))
	write(txt, file=out_file, append=TRUE)

	write.table(result, file=out_file, row.names=FALSE, col.names=FALSE, append=TRUE)
}

sfreemap.test.sfreemap <- function(tree, Q, nsim_values, hist, sample_freq) {

	result <- create_result_matrix(nsim_values)

	for (nsim in nsim_values) {
		t_start <- proc.time()
		sm <- sfreemap.map(tree
							, sample_freq=sample_freq
			  			 	, Q=Q
			  			 	, n_simulations=nsim
						 	, rewards=rep(1,ncol(tree$node.states)))
		t_elapsed <- proc.time() - t_start

		mean <- sfreemap.describe(sm)

		diff <- sfreemap.diff(hist, mean)
		data <- c(nsim, t_elapsed[3], diff$lmt, diff$emr, mean$lmt, mean$emr[1], mean$emr[2])

		result[as.character(nsim),] <- data
	}
	return(result)

}

sfreemap.test.simmap <- function(tree, Q, pi, model, nsim_values, hist, sample_freq) {

	result <- create_result_matrix(nsim_values)

	for (nsim in nsim_values) {

		t_start <- proc.time()
		mtrees <- make.simmap(tree, tree$states, Q=Q, pi=pi, model=model, nsim=nsim, samplefreq=sample_freq, message=FALSE)
		t_elapsed <- proc.time() - t_start

		mean <- simmap.mean(mtrees)

		diff <- sfreemap.diff(hist, mean)

		data <- c(nsim, t_elapsed[3], diff$lmt, diff$emr, mean$lmt, mean$emr[1], mean$emr[2])

		result[as.character(nsim),] <- data
	}

	return (result)
}

create_result_matrix <- function(nsim_values) {

	# TODO: only works for two states, should work for any number
	metric_values <- c("nsim", "time", "diff_lmt", "diff_emr", "transitions", "time_in_a", "time_in_b")
	result <- matrix(0, nrow=length(nsim_values)
			  , ncol=length(metric_values)
			  , dimnames=list(nsim_values, metric_values))

	return (result)

}
