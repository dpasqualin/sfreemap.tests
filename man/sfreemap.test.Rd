\name{sfreemap.test}
\alias{sfreemap.test}
\title{Simulation free stochastic character maps on a phylogenetic tree}
\usage{
sfreemap.test(species=100, nsim_values=NULL, Q=NULL, QS=NULL
                          , model='ER', pi="estimated", message=FALSE
                          , test_sfreemap=TRUE, test_simmap=TRUE
                          , save_to_file=TRUE
                          , ...)
}
\arguments{
    \item{species}{the number of species in the simulated phylogeny}
    \item{nsim_values}{a vector where each element represents a number of simulations}
    \item{Q}{the transition rate matrix}
    \item{QS}{the transition rate matrix used for the tree simulation}
    \item{model}{model}
    \item{pi}{the priors}
    \item{message}{TRUE if you want to see the progress on the screen}
    \item{test_sfreemap}{should the test program run for sfreemap?}
    \item{test_simmap}{should the test program run for simmap?}
    \item{save_to_file}{Should the result be written to a file?}
    \item{...}{more}
}
\description{
    Description
}
\details{
    Details
}
\value{
    Returns
}

\references{
    References
}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    See Also
}
\keyword{phylogenetics}
