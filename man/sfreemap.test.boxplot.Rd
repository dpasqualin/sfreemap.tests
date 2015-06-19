\name{sfreemap.test.boxplot}
\alias{sfreemap.test.boxplot}
\title{Simulation free stochastic character maps on a phylogenetic tree}
\usage{
sfreemap.test.boxplot(species=100
							   , Q='mcmc'
							   , pi='equal'
							   , model='ER'
							   , trees=NULL
                               , n_trees=1
							   , nsim=50
							   , n_tests=5
							   , sample_freq=100
							   , save_to_file=TRUE
			  				   , ...)
}
\arguments{
    \item{species}{the number of species in the simulated phylogeny}
	\item{Q}{The transition matrix}
	\item{pi}{whatever}
	\item{model}{The evolution model}
	\item{trees}{The trees to test on. Leave it NULL to test with a single random tree}
	\item{n_trees}{When trees=NULL you can choose how many random trees to
    generate}
    \item{nsim}{number of trees generated}
    \item{n_tests}{number of tests, used to calculate the mean per sample}
    \item{sample_freq}{the frequency of samples taken form the markov chain}
    \item{save_to_file}{whether to save the results into a file}
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