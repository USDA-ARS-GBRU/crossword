% File library/man/select.Rd
\name{select}
\alias{select}
\title{selecte}
\description{
The function selects individuals from a population.
}
\usage{
select(pop,N,level,gen,method,pheno,percentage,im_type)
}
\arguments{
\item{pop}{Required: An  input haplotype object}
\item{N}{Required: Number of individuals to be selected}
\item{level}{Required: the level of selection}
\item{method}{Required: a method of selection: random, top, or bottom}
\item{pheno}{Required: a phenotyping data, it is required with top or bottom methods}
\item{percentage}{optional: TRUE if N is percentage of a population}
\item{im_type}{Optional: the type of images, default is svg}
\item{file}{Optional: the output file path, default is current directory}
}
\value{
A haplotype object of the selected individuals
}
\examples{
haplotype_sub1 = select(haplotype,1,0,,method="random",phenotypes)}
