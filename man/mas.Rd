% File library/man/mas.Rd
\name{mas}
\alias{mas}
\title{mas}
\description{
The function selects a subset of population based on marker assistant selection supplied as input.}
\usage{
mas(pop,marker,level,parental_genotypes)
}
\arguments{
\item{pop}{Required: a haplotype object}
\item{marker}{Required: a list of marker sets and minor allel frequencies}
\item{level}{Required: a level for selection; 0, 1 or 2}
\item{parental_genotypes}{Required: a parental_genotypes object}}
\value{
A haplotype object.
}
\examples{
marker = list(c("Aradu.A04_106615237","C",0.5),c("Aradu.A03_7819678","T",0.5))
hap1 = mas(pop,marker,2,parental_genotypes)}
