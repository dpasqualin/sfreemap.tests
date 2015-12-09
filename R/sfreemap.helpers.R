summary_speed_up <- function(t1, t2) {
    ret <- list(
        'mean' = mean(apply(cbind(t1,t2), 1, function(x) x[1]/x[2]))
        , 'lastf1' = tail(t1, 1)
        , 'lastf2' = tail(t2, 1)
        , 'max' = tail(t1, 1) / tail(t2, 1)
    )
    return(ret)
}

simulation.data <- function(tree) {
    # the order of columns created by the sim.history might be different
    # from what simmap and sfreemap expect.

    if (inherits(tree, 'phylo')) {
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
    if (inherits(mtrees, 'multiPhylo')) {
        mean_emr <- rowMeans(sapply(mtrees, function(x) colSums(x$mapped.edge)))
        mean_lmt <- colMeans((countSimmap(mtrees,message=FALSE)))[1]
    } else {
        mean_emr <- colSums(mtrees$mapped.edge)
        mean_lmt <- countSimmap(mtrees)$N
    }
    return (list(lmt=mean_lmt, emr=mean_emr))
}

# This function is suppose to test whether more simulations means more accuracy
# on simmap. The results are a bit odd, maybe I'm not thinking right about
# how to compute this. Need revision..
simmap.calc_simulations_evolution <- function(trees, plot=FALSE) {
    data <- matrix(NA, nrow=length(trees), ncol=2)
    tmean <- simmap.mean(trees[[1]])$emr
    data[1,] <- c(1, 0)
    for (i in 2:length(trees)) {
        emr <- colSums(trees[[i]]$mapped.edge)
        diff <- sum(abs(emr - tmean))
        data[i,] <- c(i, diff)
        tmean <- simmap.mean(trees[1:i])$emr
    }
    data <- data.frame(data)
    colnames(data) <- c('iteration', 'diff')

    if (isTRUE(plot)) {
        breaks <- data$iteration
        p <- ggplot(data, aes(x=iteration, y=diff)) +
                geom_line() +
                scale_x_continuous(breaks=breaks) +
                theme_bw(base_size=26) +
                xlab('Tree number') +
                ylab('Distance to mean') +
                theme(axis.title.y=element_text(vjust=1.8))
        print(p)
    }

    return(data)
}

# return t1 only with tips that are in t2 too
# optionally reroot at node 'reroot'
sfreemap.pruning <- function(t1, t2, reroot=NULL) {
    tips_to_remove <- t1$tip.label[!t1$tip.label %in% t2$tip.label]
    t <- drop.tip(t1, tips_to_remove)
    if (!is.null(reroot)) {
        if (reroot %in% t$tip.label) {
            t <- root(t, reroot)
        } else {
            msg <- paste('trying to root tree but tip', reroot, 'doesn\'t exist')
            stop(msg)
        }
    }
    return(t)
}

# This function calculates the "mean tree", in other words, the mean value
# for dwelling times on states of a multiPhylo object.
sfreemap.reduce <- function(trees, type='mean') {
    reduced_trees <- list()

    start_idx <- 1
    end_idx <- 0
    cont <- 1

    keep_going <- function(t, start, end) {
        return (length(t) > end_idx && all.equal.phylo(t[[start]], t[[end]]))
    }

    while (TRUE) {

        while (keep_going(trees, start_idx, end_idx+1)) {
            end_idx <- end_idx + 1
        }

        range <- start_idx:end_idx
        current_set <- trees[range]
        base_tree <- current_set[[1]]

        if (length(range) > 1) {
            states <- colnames(base_tree$mapped.edge)
            mapped.edge <- lapply(current_set, function(x) x$mapped.edge)

            if (type == 'mean') {
                base_tree$mapped.edge <- Reduce('+', mapped.edge) / length(mapped.edge)
            } else if (type == 'median') {
                reduced <- Reduce(cbind, mapped.edge)
                for (state in states) {
                    tmp <- reduced[,colnames(reduced)==state]
                    base_tree$mapped.edge[,state] <- apply(tmp, 1, median)
                }
            } else {
                stop('unrecognized type, we only know mean and median')
            }
        }

        reduced_trees[[cont]] <- base_tree

        if (length(trees) == end_idx) {
            break
        } else {
            start_idx <- end_idx + 1
            cont <- cont + 1
        }
    }

    if (length(reduced_trees) > 1) {
        class(reduced_trees) <- 'multiPhylo'
    } else {
        reduced_trees <- reduced_trees[[1]]
        class(reduced_trees) <- 'phylo'
    }

    return(reduced_trees)
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

	write.table(result, file=out_file, row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
}

create_result_matrix <- function(n) {

	# TODO: only works for two states, should work for any number
	metric_values <- c("n_trees", "n_species", "q_size", "time", "nsim", "mode")
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

calc_time <- function(trees, parallel, prog, n_tests, n_sim, fixed_q, remove_outliers=TRUE) {

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
        if (prog == 'sfreemapc') {
            Q <- NULL
        } else {
            Q <- 'empirical'
        }
    }

    values <- rep(0, n_tests)

    for (i in 1:n_tests) {
        if (prog == 'sfreemap') {
            t <- doit(sfreemap::sfreemap.map(trees, states, Q=Q, parallel=parallel))
        } else if (prog == 'sfreemapc') {
            t <- doit(sfreemapc::sfreemap.map(trees, states, Q=Q, method='empirical', type='standard', parallel=parallel))
        } else if (prog == 'simmap') {
            t <- doit(make.simmap(trees, states, Q=Q, nsim=n_sim, message=FALSE))
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

    add_qs <- function(tree, Q) {
        tree[['Q']] <- Q
        return (tree)
    }

    if (n_trees > 1) {
        trees <- lapply(trees, add_qs, Q=QS)
        class(trees) <- 'multiPhylo'
    } else {
        trees[['Q']] <- QS
        class(trees) <- 'phylo'
    }

    return(trees)
}
