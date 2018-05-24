% File library/man/get_haplotype_matrix.Rd
\name{get_haplotype_matrix}
\alias{get_haplotype_matrix}
\title{get_haplotype_matrix}
\description{
The function reformats the individual haplotypes into a unifrom matrix
}
\usage{
get_haplotype_matrix(haplotype)
}
\arguments{
\item{haplotype}{An individual haplotype}}
\value{
The reformated haplotype in a matrix form of three columns; phyical_location, parent1 and parent2.
}
\examples{
hap1 = get_haplotype_matrix(hap1)}