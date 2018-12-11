% File library/man/combine_populations.Rd
\name{combine_populations}
\alias{combine_populations}
\title{combine_populations}
\description{
The function combines multiple haplotypes' objects into a single haplotype object \cr
It accept list of populations and merge them in single population.}
\usage{
combine_populations(pops)
}
\arguments{
\item{pops}{The first haplotype object.}}
\value{
A combined haplotype object.
}
\examples{

haplo = combine_populations(list(pop1,pop2))}
