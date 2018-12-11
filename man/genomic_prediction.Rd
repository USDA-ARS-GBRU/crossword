% File library/man/genomic_prediction.Rd
\name{genomic_prediction}
\alias{genomic_prediction}
\title{genomic_prediction}
\description{
The function applies genomic selection for the supplied prediction genotype using rrBLUP by supplied a training population ,haplotype and genotype.}
\usage{
genomic_prediction(train_geno,train_pheno,predict_geno_pop,method)
}
\arguments{
\item{train_geno}{Required: a training genotype object}
\item{train_pheno}{Required: a training phenotype object}
\item{predict_geno_pop}{Required: a testing genotype object for rrBLUP or a testing population for parental_contribution method}
\item{method}{Optional: a method of genomic prediction, rrBLUP or parental_contribution, the default is rrBLUP}}
\value{
a two-column tab delimited phenotyping object; id and value
}
\examples{
train_geno = get_genotypes(parental_genotypes,pop2) \cr
train_pheno = pheno_rand(haplo_genotypes,qtn,h2,biased_selection,gen2phy,dominant) \cr
predict_geno = get_genotypes(parental_genotypes,pop3) \cr
phen1 = genomic_prediction(train_geno,train_pheno,predict_geno)}
