% File library/man/parental_contribution.Rd
\name{parental_contribution}
\alias{parental_contribution}
\title{parental_contribution}
\description{
The function estimates phentypes by parental_contribution genomic selection method
}
\usage{
parental_contribution(pop,parental_genotypes,parental_phenotypes)
}
\arguments{
\item{pop}{Required: a population to be estimated}
\item{parental_genotypes}{Required: for training genotyping}
\item{parental_phenotypes}{Required: for training phenotyping}}
\value{
a two-column tab delimted phenotyping object; id and value
}
\examples{
pheno5 = parental_contribution(pop3,parental_genotypes,parental_phenotypes)
}
