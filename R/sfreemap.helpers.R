simulation.data <- function(tree) {
    # the order of columns created by the sim.history might be different
    # from what simmap and sfreemap expect.

    if (class(tree) == 'phylo') {
        fix_order <- tree$mapped.edge[,order(colnames(tree$mapped.edge))]
        lmt <- countSimmap(tree)$N
        emr <- apply(fix_order, 2, sum)
    } else {
        lmt <- mean(countSimmap(tree)$Tr[,1])
        emr <- t(sapply(tree, function(x) colSums(x$mapped.edge)))
        emr <- colMeans(emr)
        emr <- emr[order(names(emr))]
    }
    return (list(lmt=lmt, emr=emr))
}

sfreemap.diff <- function(a, b) {
    lmt <- abs(a$lmt - b$lmt)
    emr <- sum(abs(a$emr - b$emr))
    return (list(lmt=lmt, emr=emr))
}

simmap.mean <- function(mtrees) {
    if (class(mtrees) == 'multiPhylo') {
        mean_emr <- apply(sapply(mtrees, function(x) { apply(x$mapped.edge,2,sum)}), 1, mean)
        mean_lmt <- colMeans((countSimmap(mtrees,message=FALSE)))[1]
    } else {
        mean_emr <- apply(mtrees$mapped.edge,2,sum)
        mean_lmt <- countSimmap(mtrees)$N
    }
    return (list(lmt=mean_lmt, emr=mean_emr))
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

write_to_file <- function(out_file, result, tree=NULL, out_dir=NULL
                            , hist=NULL) {

    if (!is.null(tree)) {
    	save_tree_file(out_dir, tree)
    }

    if (!is.null(hist)) {
    	txt <- paste('# sim.history ', 0, hist$lmt
                        , paste(hist$emr, collapse=' '))
    	write(txt, file=out_file)
    }

	txt <- paste('#',paste(colnames(result), collapse=','))
	write(txt, file=out_file, append=TRUE)

	write.table(result, file=out_file, row.names=FALSE, col.names=FALSE, append=TRUE)
}

create_result_matrix <- function(n) {

	# TODO: only works for two states, should work for any number
	metric_values <- c("n_trees", "n_species", "q_size", "time")
	result <- matrix(0, nrow=n
			  , ncol=length(metric_values)
			  , dimnames=list(1:n, metric_values))

	return (result)
}

calc_time <- function(trees, parallel, prog, n_tests) {
    if (class(trees)=='phylo') {
        states <- trees$states
    } else {
        states <- trees[[1]]$states
    }
    t_start <- proc.time()
    for (i in 1:n_tests) {
        if (prog == 'sfreemap') {
            invisible(sfreemap.map(trees
                                    , states
                                    , Q='empirical'
                                    , parallel=parallel))
        } else if (prog == 'simmap') {
            invisible(make.simmap(trees, states, Q='empirical'))
        } else {
            stop('prog should be equal to sfreemap or simmap')
        }
    }
    t_elapsed <- (proc.time() - t_start)[3]/n_tests
    return(t_elapsed)
}

create_trees <- function(n_trees, n_species, q_size) {
    # Create Q
    QS <- matrix(1, nrow=q_size, ncol=q_size)
    diag(QS) <- -sum(QS[1,]) + 1
    rownames(QS) <- colnames(QS) <- 1:nrow(QS)

    # Create topologies
    topologies <- pbtree(n=n_species, nsim=n_trees, scale=1)
    if (n_trees > 1) {
        trees <- lapply(topologies, sim.history, Q=QS, message=FALSE)
        class(trees) <- 'multiPhylo'
    } else {
        trees <- sim.history(topologies, QS, message=FALSE)
    }
    return(trees)
}
