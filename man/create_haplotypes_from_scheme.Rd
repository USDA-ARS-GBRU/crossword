% File library/man/create_haplotypes_from_scheme.Rd
\name{create_haplotypes_from_scheme}
\alias{create_haplotypes_from_scheme}
\title{create_haplotypes_from_scheme}
\description{
The function simulates haplotypes from a scheme \cr
It produces a haplotype object following the individuals and generations of a scheme file.}
\usage{
create_haplotypes_from_scheme(scheme,chr_length,pop)
}
\arguments{
\item{scheme}{Required: an input scheme}
\item{chr_length}{Required: a chr_length file}
\item{pop}{Optional: an input population}}
\value{
haplotypes
}
\examples{
haplotype2 = create_haplotypes_from_scheme(scheme,"chr_length.txt",pop2)
}
