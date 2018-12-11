% File library/man/haplotypes_out.Rd
\name{haplotypes_out}
\alias{haplotypes_out}
\title{haplotypes_out}
\description{
The function creates a readable matrix for input haplotypes and/or store it for further analysis.}
\usage{
haplotypes_out(pop,output,parental_genotypes)
}
\arguments{
\item{pop}{Required: an input haplotype object}
\item{output}{Optional: Output prefix for the csv and haplo files}
\item{parental_genotypes}{Required: parental genotyping object}}
\examples{
haplotypes_out(haplotypes,"output",parental_genotypes)}
