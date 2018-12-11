% File library/man/get_values.Rd
\name{get_values}
\alias{get_values}
\title{get_values}
\description{
The function calculates the phenotypic values based on th supplied input genotype object and QTN effects.}
\usage{
get_values(geno,qtn_effect,dominant)
}
\arguments{
\item{geno}{Required: genotype object}
\item{qtn_effect}{Required: qtn_effect object}
\item{dominant}{Optional: a logical object for dominance, default is FALSE}}
\value{
the phenotypic values
}
\examples{
gff = "crossword_processing/peanut.gff"
chr_stat = "crossword_processing/chr_siz.txt"
chr_length = "crossword_processing/chr_len.txt"
gen2phy = physical2genomic(gff,chr_stat,chr_length,200000)
homo = TRUE
input_pheno = "crossword_processing/pheno.txt"
parental_genotypes = get_parental_genotypes(input,gen2phy,homo)
biased_selection = TRUE
min_qtn_freq = 0
dominant=FALSE
effect_distribution = "equall"
pop3 = create_families(pop=pop2,S=15,id="pop3",chr_length=chr_length)
ids = pop3[[1]][pop3[[1]]$gen==pop3[[1]][nrow(pop3[[1]]),2],1]
geno = select_haplotype(pop3,ids)
qtn_effect = get_values(qtn = 5,gen2phy=gen2phy,biased_selection=biased_selection,parental_genotypes=parental_genotypes,min_qtn_freq=min_qtn_freq,dominant=min_qtn_freq,effect_distribution=effect_distribution)
values = get_values(geno,qtn_effect,dominant)}
