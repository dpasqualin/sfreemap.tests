
create_plots <- function(lang='pt_BR') {

    # this is by far not the best way to do it, but I feel lazy right now
    if (lang == 'pt_BR') {
        states <- "Número de estados"
        taxa <- "Número de taxa"
        trees <- "Número de árvores"
        cores <- "Número de núcleos de processamento"
        serial <- "Serial"
        parallel <- "Paralelo"
    } else {
        states <- "Number of states"
        taxa <- "Number of taxa"
        trees <- "Number of trees"
        cores <- "Number of cores"
        serial <- "Serial"
        parallel <- "Parallel"
    }


    x <- list(
        files=c('sfreemapc.states.txt', 'sfreemap.states.txt', 'simmap.states.txt'),
        types=c('state', 'state', 'state'),
        time=c('time_to_estimate', 'time_to_estimate', 'time_to_estimate'),
        legend=c('SFREEMAP-C', 'SFREEMAP-R', 'SIMMAP')
    )

    plot_comparison_for_q(x, states, trans='log10', output='/tmp/estimate_q_states.png')
    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.taxa.txt', 'sfreemap.taxa.txt', 'simmap.taxa.txt'),
        types=c('taxa', 'taxa', 'taxa'),
        time=c('time_to_estimate', 'time_to_estimate', 'time_to_estimate'),
        legend=c('SFREEMAP-C', 'SFREEMAP-R', 'SIMMAP')
    )

    plot_comparison_for_q(x, taxa, output='/tmp/estimate_q_taxa.png')

    # ----------------------------------------------------------
    #x <- list(
    #    files=c('sfreemapc.states.txt', 'sfreemap.states.txt', 'sfreemapc.states.txt', 'sfreemap.states.txt'),
    #    types=c('state', 'state', 'state', 'state'),
    #    time=c('time_to_map', 'time_to_map', 'time_to_estimate', 'time_to_estimate'),
    #    legend=c('SFREEMAP-C-MAP', 'SFREEMAP-R-MAP', 'SFREEMAP-C-ESTIMATE', 'SFREEMAP-R-ESTIMATE')
    #)

    #plot_comparison_for_q(x, states, output='/tmp/analyse_q_states.png')

    # ----------------------------------------------------------
    x <- list(
        files=c('sfreemapc.trees.txt', 'sfreemap.trees.txt'),
        types=c('tree', 'tree'),
        legend=c('SFREEMAP-C', 'SFREEMAP-R')
    )

    plot_comparison(x, trees, lang=lang, output='/tmp/trees_serial.png')

    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.trees.txt', 'sfreemapc.trees.txt'),
        types=c('tree', 'tree'),
        mode=c('serial', 'parallel'),
        legend=c(serial, parallel)
    )

    plot_comparison(x, trees, lang=lang, output='/tmp/trees_parallel.png')

    # ----------------------------------------------------------
    x <- list(
        files=c('sfreemapc.taxa.txt', 'sfreemap.taxa.txt'),
        types=c('taxa', 'taxa'),
        legend=c('SFREEMAP-C', 'SFREEMAP-R')
    )

    plot_comparison(x, taxa, lang=lang, output='/tmp/taxa_serial.png')

    # ----------------------------------------------------------
    x <- list(
        files=c('sfreemapc.states.txt', 'sfreemap.states.txt'),
        types=c('state', 'state'),
        legend=c('SFREEMAP-C', 'SFREEMAP-R')
    )

    plot_comparison(x, states, lang=lang, trans="log10", output='/tmp/states_serial.png')
    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.states.txt', 'simmap.states.txt', 'simmap.states.txt', 'simmap.states.txt'),
        types=c('state', 'state', 'state', 'state'),
        nsim=c(1, 1, 10, 20),
        legend=c('SFREEMAP', 'SIMMAP-1', 'SIMMAP-10', 'SIMMAP-20')
    )

    plot_comparison(x, states, lang=lang, trans="log10", output='/tmp/states_simmap.png')
    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.taxa.txt', 'simmap.taxa.txt', 'simmap.taxa.txt', 'simmap.taxa.txt'),
        types=c('taxa', 'taxa', 'taxa', 'taxa'),
        nsim=c(1, 1, 10, 20),
        legend=c('SFREEMAP', 'SIMMAP-1', 'SIMMAP-10', 'SIMMAP-20')
    )

    plot_comparison(x, taxa, lang=lang, output='/tmp/taxa_simmap.png')
    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.trees.txt', 'simmap.trees.txt', 'simmap.trees.txt', 'simmap.trees.txt'),
        types=c('tree', 'tree', 'tree', 'tree'),
        nsim=c(1, 1, 10, 20),
        legend=c('SFREEMAP', 'SIMMAP-1', 'SIMMAP-10', 'SIMMAP-20')
    )

    plot_comparison(x, trees, lang=lang, output='/tmp/trees_simmap.png')
    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.states.omp.txt', 'sfreemapc.states.omp.txt'),
        types=c('omp', 'omp'),
        q=c('fixed', 'estimated'),
        legend=c('Mapeamento', 'Mapeamento mais matriz Q')
    )

    plot_comparison(x, cores, lang=lang, output='/tmp/states_omp.png')
    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.states.omp.txt', 'sfreemapc.states.omp.txt'),
        types=c('omp', 'omp'),
        q=c('fixed', 'estimated'),
        legend=c('Mapeamento', 'Mapeamento mais matriz Q')
    )

    plot_speed_up(x, output='/tmp/states_speedup_omp.png')

    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.trees.mccores.txt'),
        types=c('cores'),
        mode=c('parallel'),
        q=c('estimated'),
        legend=c('Sfreemap Paralelo')
    )

    plot_speed_up(x, limit=17, output='/tmp/trees_parallel_speedup.png')

    return(NULL)
}
