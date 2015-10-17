# TODO: just starting this function, not sure if it sue be here or in
# sfreemap.tests package.
# Based on the plot found in http://www.inside-r.org/r-doc/graphics/pairs
sfreemap.plot_correlation <- function(sfreemap_empirical, sfreemap_mcmc, simmap_empirical, simmap_mcmc, state='1') {
    # panel.smooth function is built in.
    # panel.cor puts correlation in upper panels, size proportional to correlation
    # from: http://www.r-bloggers.com/scatterplot-matrices-in-r/
    panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits=digits)[1]
        txt <- paste(prefix, txt, sep="")
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
    }

    # add histogram on the diagonals
    # from: https://stat.ethz.ch/R-manual/R-devel/library/graphics/html/pairs.html
    panel.hist <- function(x, ...) {
        print(x)
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5))
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$counts
        y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
    }

    # build data.frame
    data <- cbind(sfreemap_empirical$mapped.edge[,state]
                    , sfreemap_mcmc$mapped.edge[,state]
                    , simmap_empirical$mapped.edge[,state]
                    , simmap_mcmc$mapped.edge[,state])
    rownames(data) <- NULL
    colnames(data) <- c('sfreemap_empirical', 'sfreemap_mcmc'
                        , 'simmap_empirical', 'simmap_mcmc')
    data <- data.frame(data, stringsAsFactors=FALSE)

    # add color information to it
    nodes <- as.numeric(sapply(strsplit(rownames(a), ','), function(x) x[2]))
    nodes <- tmp <- sfreemap_empirical$tip.label[nodes]
    nodes[is.na(tmp)] <- 'grey'
    nodes[!is.na(tmp)] <- 'green'
    data$color <- nodes

    # package car has a scatterplotMatrix function that can print things like
    # an ellipse showing 95% confidence level, but I'm not sure it this fits
    # here. https://cran.r-project.org/web/packages/car/car.pdf
    pairs(~ sfreemap_empirical + sfreemap_mcmc + simmap_empirical + simmap_mcmc
        , data=data
        , pch = 21
        , bg = data$color
        , upper.panel=panel.cor
        , diag.panel=panel.hist
    )

    return(data)
}
