% File library/man/create_population.Rd
\name{create_population}
\alias{create_population}
\title{create_population}
\description{
creating a population from parents
}
\usage{
create_population(P,id,chr_length,parental_genotypes)
}
\arguments{
\item{P}{Required: a list of parents to be used}
\item{chr_length}{Required: a two-column file of chromosomes' IDs and length in CM}
\item{parental_genotypes}{Required: parental_genotypes object}}
\examples{
 pop1 = create_population(P= c("tr","tg","grg"),chr_length,parental_genotypes)}
