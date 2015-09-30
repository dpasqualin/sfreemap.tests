\name{sfreemap.test.boxplot}
\alias{sfreemap.test.boxplot}
\title{Simulation free stochastic character maps on a phylogenetic tree}
\usage{
sfreemap.test.boxplot(species=100
                        , Q_simmap='mcmc'
                        , Q_sfreemap='empirical'
					    , trees=NULL
                       	, n_trees=1
					    , nsim=50
				        , n_tests=5
					    , sample_freq=100
					    , save_to_file=TRUE
	  				    , ...)
}
\arguments{
    \item{species}{The number of species in the simulated phylogeny}
    \item{Q_mcmc}{The transition rate matrix for simmap}
	\item{Q_sfreemap}{The transition rate matrix for sfreemap}
	\item{trees}{The trees to test on. Leave it NULL to test with random trees}
	\item{n_trees}{When trees=NULL you can choose how many random trees to
    generate}
    \item{nsim}{...}
    \item{n_tests}{number of tests to run, used to calculate the average}
    \item{sample_freq}{the frequency of samples taken form the markov chain}
    \item{save_to_file}{a boolean indicating whether to save the results into a file}
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
