% File library/man/genotypes_out.Rd
\name{genotypes_out}
\alias{genotypes_out}
\title{genotypes_out}
\description{
The function exports genotypes into hapmap format and SNP genomic locations in a three column format (chr, pos, and loci).
}
\usage{
genotypes_out(haplo_genotypes,output,pop,parental_genotypes,A_parent, B_parent)
}
\arguments{
\item{haplo_genotypes}{Required: an input genotypes object, the output of "get_genotypes" function}
\item{output}{Required: a prefix of the output files}
\item{level}{Optional: the level of outputs}
\item{pop}{Optional: the input population, required if level is assigned }
\item{parental_genotypes}{Optional: the parental genotype object, is required if the user needs output in map format of biparental}
\item{A_parent}{Optional: the parent A in mapping of bipaerntal}
\item{B_parent}{Optional: the parent B in mapping of bipaerntal}}
\examples{

genotypes_out(haplo_genotypes,"haplo_genotypes.csv")}
