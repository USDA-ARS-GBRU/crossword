% File library/man/get_marker_intersection.Rd
\name{get_marker_intersection}
\alias{get_marker_intersection}
\title{get_marker_intersection}
\description{
The function creates marker sub sets from genotpe object, three options are available:  \cr
1. an input file of a set of markers \cr
2. a phenotype object \cr
3. another genotype object
}
\usage{
get_marker_intersection(genotypes,input_geno_file,pheno,geno)
}
\arguments{
\item{genotypes}{Required: input genotype object}
\item{input_geno_file}{Optional: an input genotype file}
\item{pheno}{Optional: a phenotype object}
\item{geno}{Optional: a genotype object}}
\value{
a genotype object of two datasets; genotypes and gen2phy
}
\examples{
get_marker_intersection(genotypes=genotypes,input_geno_file = input_geno_file)}
