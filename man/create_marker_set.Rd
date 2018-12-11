% File library/man/create_marker_set.Rd
\name{create_marker_set}
\alias{create_marker_set}
\title{create_marker_set}
\description{
The function creates marker sub sets from genotype object of N number of random marker with a certain allele frequency cutoff.}
\usage{
create_marker_set(genotypes,N,MAF,input_geno_file,pheno,geno,biased_selection,gen2phy,parental_genotypes)
}
\arguments{
\item{genotypes}{Required: input genotype object}
\item{N}{Optional: number of markers to be selected randomly}
\item{MAF}{Optional: minor allels frequency cutoff}
\item{no_qtn}{Optional: a list of qtn to be excluded from the selected markers}
\item{biased_selection}{optional: if TRUE, the randomly selected QTNs are biased to gene denesity regions}
\item{gen2phy}{optinal: gen2phy object to calculate biased_selection}
\item{parental_genotypes}{Required: parental_genotypes object}}
\value{
a genotype object of two datasets; genotypes and gen2phy
}
\examples{
ms1 = create_marker_set(genotypes=genotypes,N=6,MAF=0.2)}
