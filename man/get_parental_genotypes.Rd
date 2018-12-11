% File library/man/get_parental_genotypes.Rd
\name{get_parental_genotypes}
\alias{get_parental_genotypes}
\title{get_parental_genotypes}
\description{
The function collects the SNPs and their genotypes/gen2phy information for parental genotypes.}
\usage{
get_parental_genotypes(input,gen2phy,keep_homo)
}
\arguments{
\item{input}{Required: a input vcf or hapmap file.}
\item{gen2phy}{Required: gen2phy object, the output of "physical2genomic" function.}
\item{keep_homo}{Optional: if TRUE only homozygous paretnal loci will be keept, otherwise the output will have homozygous and heterozygous parental loci}}
\value{
A structured output of: \cr
genotypes: the parental genotypes of SNPs. \cr
gen2phy: the genetic loci of the SNPs
}
\examples{
parental_genotypes = get_parental_genotypes('three_genotypes_snps.vcf',gen2phy,TRUE)}
