% File library/man/create_families.Rd
\name{create_families}
\alias{create_families}
\title{create_families}
\description{
The function create segregant families from a population \cr
It takes a cross, create_population or advance output object as input and creates lines of supplied replicates from every single individual in the last generation in the input object.}
\usage{
create_families(pop,S,id,chr_length,is_clone)
}
\arguments{
\item{pop}{Required: a cross population}
\item{S}{Required: Number of individuals within families}
\item{id}{Required: the new population ID}
\item{chr_length}{Required: a file contains lengths of chromosomes}
\item{is_clone}{Optional: TRUE if the families are clones, default is FALSE}}
\value{
A haplotype object of N individual with lowest phenotypic breeding value
}
\examples{
pop3 = create_families(pop=pop1,S=15,id="pop3",chr_length=chr_length)}
