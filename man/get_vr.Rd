% File library/man/get_vr.Rd
\name{get_vr}
\alias{get_vr}
\title{get_vr}
\description{
The function calculate Vr for a certain supplied heritability.}
\usage{
get_vr(qtn_effect,h2,parental_genotypes,dominant,heritability_mode)
}
\arguments{
\item{qtn_effect}{Required: QTN effect matrix}
\item{h2}{optinal: heritability, default is 0.5}
\item{parental_genotypes}{Required: parental_genotypes object}
\item{dominant}{Optional: a logical object for dominance, default is FALSE}
\item{heritability_mode}{Required: heritability_mode is "absolute", "relative", or "pick_parent"}}
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
qtn_effect = random_qtn_assign(qtn = 5,gen2phy=gen2phy,biased_selection=biased_selection,parental_genotypes=parental_genotypes,min_qtn_freq=min_qtn_freq,dominant=min_qtn_freq,effect_distribution=effect_distribution)
h2 = 0.6
tbv_only = FALSE
vr1 = get_vr(qtn_effect = qtn_effect,h2 = 0.6 ,parental_genotypes=parental_genotypes,dominant=TRUE,heritability_mode="absolute")
}
