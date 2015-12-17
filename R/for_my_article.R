
create_plots <- function() {

    # ----------------------------------------------------------
    x <- list(
        files=c('sfreemapc.states.txt', 'sfreemap.states.txt', 'sfreemapc.states.txt', 'sfreemap.states.txt'),
        types=c('state', 'state', 'state', 'state'),
        time=c('time_to_map', 'time_to_map', 'time_to_estimate', 'time_to_estimate'),
        legend=c('SFREEMAP-C-MAP', 'SFREEMAP-R-MAP', 'SFREEMAP-C-ESTIMATE', 'SFREEMAP-R-ESTIMATE')
    )

    plot_comparison_for_q(x, "Número de estados", output='/tmp/analyse_q_states.png')

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
        files=c('sfreemapc.states.txt', 'simmap.states.txt', 'simmap.states.txt', 'simmap.states.txt'),
        types=c('state', 'state', 'state', 'state'),
        nsim=c(1, 1, 10, 20),
        legend=c('SFREEMAP', 'SIMMAP-1', 'SIMMAP-10', 'SIMMAP-20')
    )

    plot_comparison(x, "Número de estados do caráter", output='/tmp/states_simmap.png')
    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.taxa.txt', 'simmap.taxa.txt', 'simmap.taxa.txt', 'simmap.taxa.txt'),
        types=c('taxa', 'taxa', 'taxa', 'taxa'),
        nsim=c(1, 1, 10, 20),
        legend=c('SFREEMAP', 'SIMMAP-1', 'SIMMAP-10', 'SIMMAP-20')
    )

    plot_comparison(x, "Número de taxa", output='/tmp/taxa_simmap.png')
    # ----------------------------------------------------------

    x <- list(
        files=c('sfreemapc.trees.txt', 'simmap.trees.txt', 'simmap.trees.txt', 'simmap.trees.txt'),
        types=c('tree', 'tree', 'tree', 'tree'),
        nsim=c(1, 1, 10, 20),
        legend=c('SFREEMAP', 'SIMMAP-1', 'SIMMAP-10', 'SIMMAP-20')
    )

    plot_comparison(x, "Número de árvores", output='/tmp/trees_simmap.png')
    # ----------------------------------------------------------

    return(NULL)
}
