% File library/man/cross.Rd
\name{cross}
\alias{cross}
\title{cross}
\description{The function is used to create haplotype information based on the input \cr
it uses "create_scheme_from_single_population" and "create_scheme" local functionis to create the scheme \cr
and "create_parent" and "cross" functions of 'simcross' package \cr
The function creates scheme from the provided parental information: \cr
1. single parent1 X single parent2 \cr
2. single parent1 X list of parents2 \cr
3. list of parents1 X single parent2 \cr
 combinations of list of parents1 X list of parents2 is not allowed}
\usage{
cross(P1,P2,N,id,chr_length,parental_genotypes,isdh)
}
\arguments{
\item{P1}{Required: parent/s 1, a name or a list of names of the parents 1}
\item{P2}{Optional: parent/s 2, a name or a list of names of the parents 2 \cr
If it is missing, the multiple cross will be applied for all parents listed in parent 1 which should be > two parents}
\item{N}{Required: number of individuals in each repeated generation}
\item{id}{Required: ID of the current breeding procedure}
\item{chr_length}{Required: a two-column file of chromosomes' IDs and length in CM}
\item{parental_genotypes}{Required: a parental_genotypes object}
\item{isdh}{Optional: is double habloid ?}}
\value{
The output is a structured object of:
1. scheme: contains the pedigree information of all haplotypes in a data frame of four columns: \cr
	1.haplotype IDs \cr
	2.generation \cr
	3.parent 1 \cr
	4.parent 2 \cr
2. haplotype: contains the haplotype genetic architect for all haplotypes \cr
	each haplotype has sub datasets of chromosomes and each chromosome has the maternal (mat) and paternal (pat) information which is allele structure and loci
}
\examples{
chr_length = "chr_len.txt"
parental_genotypes = get_parental_genotypes(input,gen2phy,homo)
pop1 = cross("grg","tg",400,'pop1',chr_length,parental_genotypes)}
