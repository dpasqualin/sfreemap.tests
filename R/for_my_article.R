summary_speed_up <- function(f1, f2, filter='serial') {
    t1 <- read.table(f1)
    t1 <- t1[t1$V6 == filter,]$V4

    t2 <- read.table(f2)
    t2 <- t2[t2$V6 == filter,]$V4

    ret <- list(
        'mean' = mean(apply(cbind(t1,t2), 1, function(x) x[1]/x[2]))
        , 'last f1' = tail(t1, 1)
        , 'last f2' = tail(t2, 1)
        , 'speed up max' = tail(t1, 1) / tail(t2, 1)
    )

    return(ret)
}


just_do_it <- function() {

    # generate grouwth rate chart
    #x <- list(
    #    files=c('sfreemapc.taxa.txt', 'sfreemapc.states.txt',
    #            'sfreemapc.trees.txt', 'simmap.states.txt'),
    #    types=c('taxa', 'state', 'tree', 'state'),
    #    legend=c('Sfreemap-Taxa', 'Sfreemap-Estados',
    #             'Sfreemap-Árvores', 'Simmap-Estados')
    #)

    # limit=10 because we only have 10 values for states
    #plot_growth_rate(x, 10, output='/tmp/growth.png')

    # ----------------------------------------------------------
    x <- list(
        files=c('sfreemapc.trees.txt', 'sfreemap.trees.txt'),
        types=c('tree', 'tree'),
        legend=c('SFREEMAP-C', 'SFREEMAP-R')
    )

    plot_comparison(x, "Número de árvores", output='/tmp/trees_serial.png')

    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.trees.txt', 'sfreemapc.trees.txt'),
        types=c('tree', 'tree'),
        mode=c('serial', 'parallel'),
        legend=c('Serial', 'Paralelo')
    )

    plot_comparison(x, "Número de árvores", output='/tmp/trees_parallel.png')

    # ----------------------------------------------------------
    x <- list(
        files=c('sfreemapc.taxa.txt', 'sfreemap.taxa.txt'),
        types=c('taxa', 'taxa'),
        legend=c('SFREEMAP-C', 'SFREEMAP-R')
    )

    plot_comparison(x, "Número de taxa", output='/tmp/taxa_serial.png')

    # ----------------------------------------------------------
    x <- list(
        files=c('sfreemapc.states.txt', 'sfreemap.states.txt'),
        types=c('state', 'state'),
        legend=c('SFREEMAP-C', 'SFREEMAP-R')
    )

    plot_comparison(x, "Número de estados do caráter", output='/tmp/states_serial.png')
    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.states.txt', 'simmap.states.txt'),
        types=c('state', 'state'),
        legend=c('SFREEMAP', 'SIMMAP')
    )

    plot_comparison(x, "Número de estados do caráter", output='/tmp/states_simmap.png')
    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.taxa.txt', 'simmap.taxa.txt'),
        types=c('taxa', 'taxa'),
        legend=c('SFREEMAP', 'SIMMAP')
    )

    plot_comparison(x, "Número de taxa", output='/tmp/taxa_simmap.png')
    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.trees.txt', 'simmap.trees.txt'),
        types=c('tree', 'tree'),
        legend=c('SFREEMAP', 'SIMMAP')
    )

    plot_comparison(x, "Número de árvores", output='/tmp/trees_simmap.png')
    # ----------------------------------------------------------

    return(NULL)
}
