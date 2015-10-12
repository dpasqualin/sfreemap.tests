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

simmap.mean_tree <- function(trees, nsim, samplefreq=100) {
    if (is.null(trees[[1]]$maps)) {
        # calculate first and than plot
    }

    total <- length(trees)/nsim
    mean_trees <- list()
    for (i in 1:total) {
        start <- ((i-1)*nsim+1)
        end <- (i*nsim)
        range <- start : end
        t <- trees[[start]]
        # remember that t$maps will be wrong here, but it's ok because we
        # don't need it here
        mapped.edge <- lapply(trees[range], function(x) x$mapped.edge)
        t$mapped.edge <- Reduce('+', mapped.edge) / length(mapped.edge)
        mean_trees[[i]] <- t
    }
    class(mean_trees) <- 'multiPhylo'

    return(mean_trees)
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

plot_boxplot <- function(out_dir, out_file, data, y, xlabel, ylabel, line_data) {

    data <- data.frame(data)

	output <- paste(out_dir, out_file, sep='/')

	png(output, width=1024, height=768)
    p <- ggplot(data, aes_string(x='factor(generation)', y=y)) +
            theme_bw(base_size=26) +
            geom_boxplot() +
            xlab(xlabel) +
            ylab(ylabel) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            theme(axis.title.y=element_text(vjust=1.8)) +
            geom_hline(yintercept=line_data, color='red')
    print(p)
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
	metric_values <- c("n_trees", "n_species", "q_size", "time", "nsim")
	result <- matrix(0, nrow=n
			  , ncol=length(metric_values)
			  , dimnames=list(1:n, metric_values))

	return (result)
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y[!y %in% NA]
}

calc_time <- function(trees, parallel, prog, n_tests, n_sim, omp, remove_outliers=TRUE) {
    if (class(trees)=='phylo') {
        states <- trees$states
    } else {
        states <- trees[[1]]$states
    }

    if (isTRUE(omp)) {
        omp <- detectCores()
    } else {
        omp <- 1
    }

    doit <- function(expr) {
        t_start <- proc.time()
        expr
        t_end <- proc.time()
        return ((t_end-t_start)[3])
    }

    values <- rep(0,n_tests)

    for (i in 1:n_tests) {
        if (prog == 'sfreemap') {
            t <- doit(sfreemap::sfreemap.map(trees, states, Q='empirical'
                                , parallel=parallel))
        } else if (prog == 'sfreemapc') {
            t <- doit(sfreemapc::sfreemap.map(trees, states, Q='empirical'
                                , parallel=parallel, omp=omp))
        } else if (prog == 'simmap') {
            t <- doit(make.simmap(trees, states, Q='empirical'
                                  , nsim=n_sim, message=FALSE))
        } else {
            stop('valid for "prog": (simmap|sfreemap|sfreemapc)')
        }
        values[i] <- t
    }

    if (isTRUE(remove_outliers)) {
        values <- remove_outliers(values)
    }

    t_elapsed <- mean(values)

    return(t_elapsed)
}

create_trees <- function(n_trees, n_species, q_size, unique=FALSE) {
    # Create Q
    QS <- matrix(1, nrow=q_size, ncol=q_size)
    diag(QS) <- -sum(QS[1,]) + 1
    rownames(QS) <- colnames(QS) <- 1:nrow(QS)

    # Create topologies
    if (isTRUE(unique)) {
        topologies <- pbtree(n=n_species, nsim=n_trees, scale=1)
        trees <- lapply(topologies, sim.history, Q=QS, message=FALSE)
    } else {
        topology <- pbtree(n=n_species, nsim=1, scale=1)
        trees <- sim.history(topology, Q=QS, message=FALSE)
        if (n_trees > 1) {
            trees <- rep(trees, n_trees)
        }
    }

    if (n_trees > 1) {
        class(trees) <- 'multiPhylo'
    } else {
        class(trees) <- 'phylo'
    }

    return(trees)
}
