% File library/man/phenotype_out.Rd
\name{phenotype_out}
\alias{phenotype_out}
\title{phenotype_out}
\description{
The function exports genotypes into hapmap format and SNP genomic locations in a three column format (chr, pos, and loci).
}
\usage{
phenotype_out(pheno,file,im_type)
}
\arguments{
\item{pheno}{Required: an input phenotype object}
\item{file}{Required: a prefix of the output file}
\item{im_type}{Optional: the type of images, default is svg}}
\examples{

phenotype_out(pheno,"phenotype1.csv")}
