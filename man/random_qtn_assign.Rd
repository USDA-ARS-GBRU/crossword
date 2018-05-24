% File library/man/random_qtn_assign.Rd
\name{random_qtn_assign}
\alias{random_qtn_assign}
\title{random_qtn_assign}
\description{
The function select random QTNs with a certain criteria
}
\usage{
random_qtn_assign(qtn,gen2phy,biased_selection,parental_genotypes,min_qtn_freq,dominant,effect_distribution,highest_P,lowest_P,high_to_low_percentage)
}
\arguments{
\item{qtn}{Required: Number of QTNs to be selected}
\item{gen2phy}{optinal: gen2phy object to calculate biased_selection}
\item{biased_selection}{optional: if TRUE, the randomly selected QTNs are biased to gene denesity regions}
\item{parental_genotypes}{Required: parental_genotypes object}
\item{min_qtn_freq}{Optional: minimum QTN frequency, dafault is 0}
\item{dominant}{Optional: a logical object for dominance, default is FALSE}
\item{h2}{Optional: effect_distribution is "equal", "gamma", or "normal", default is "equal"}
\item{highest_P}{Optional: the highest parent in case of using pick_parent method}
\item{lowest_P}{Optional: the lowest parent in case of using pick_parent method}
\item{high_to_low_percentage}{Optional: the perecentage of opposite effect of lower parent, default 0}}
\value{
ranomly assigned QTNs and effects
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
effect_distribution = "equal"
qtn_effect = random_qtn_assign(qtn = 5,gen2phy=gen2phy,biased_selection=biased_selection,parental_genotypes=parental_genotypes,min_qtn_freq=min_qtn_freq,dominant=min_qtn_freq,effect_distribution=effect_distribution)}
