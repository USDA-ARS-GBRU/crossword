% File library/man/pick_individual.Rd
\name{pick_individual}
\alias{pick_individual}
\title{pick_individual}
\description{
The function picks an individual from a population from certain supplied level.}
\usage{
pick_individual(pop,cross,family,individual)
}
\arguments{
\item{pop}{Required: an input population}
\item{cross}{Optional: the cross level to be select from}
\item{family}{Optional: the family level to be select from}
\item{individual}{Required: an individual to be selected}}
\value{
A haplotype object of the selected individuals
}
\examples{
individual_1 = pick_individual(pop1,2,1,3)}
