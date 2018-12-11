% File library/man/get_genotypes.Rd
\name{get_genotypes}
\alias{get_genotypes}
\title{get_genotypes}
\description{
The function calculates the genotypes of certain haplotypes, supplied as an input, based on the parental genotypes.}
\usage{
get_genotypes(parental_genotypes,pop,pre_selected_markers)
}
\arguments{
\item{parental_genotypes}{Required: the parental genotypes object, the output of "get_parental_genotypes" function.}
\item{pop}{Required: a haplotypes input object, the output of "create_haplotypes" function or any other selection subset haplotype functions}
\item{pre_selected_markers}{Optional: if supplied, only these markers will be returned, it may be a file name or another genotype}
}
\value{
A structured list of a data set of haplotypes' genotypes and a data set of genomic loci
}
\examples{
pop_genotypes = get_genotypes(parental_genotypes,pop1)}

