% File library/man/phenotype.Rd
\name{phenotype}
\alias{phenotype}
\title{phenotype}
\description{
The function calculates phenotying
}
\usage{
phenotype(pop,qtn_effect,tgv_only,vr,parental_genotypes)
}
\arguments{
\item{pop}{Required: a population}
\item{qtn_effect}{Required: qtn_effect data frame}
\item{tgv_only}{Required: if TRUE, only calculated breeding value will be export, without environmental effects, default is FALSE}
\item{vr}{Required: Vr value}
\item{parental_genotypes}{Required: parental_genotypes object}}
\value{
Vr value
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
qtn_effect = random_qtn_assign(qtn = 5,gen2phy=gen2phy,biased_selection=biased_selection,parental_genotypes=parental_genotypes,min_qtn_freq=min_qtn_freq,dominant=min_qtn_freq,effect_distribution=effect_distribution)
h2 = 0.6
tgv_only = FALSE

vr1 = get_vr(qtn_effect = qtn_effect,h2 = 0.6 ,parental_genotypes=parental_genotypes,dominant=TRUE,heritability_mode="absolute")
phen1 = get_vr(pop=pop3,qtn_effect=qtn_effect,vr=vr1)
}
