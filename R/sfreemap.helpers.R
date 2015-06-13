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
