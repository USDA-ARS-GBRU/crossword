% File library/man/select_haplotype.Rd
\name{select_haplotype}
\alias{select_haplotype}
\title{select_haplotype}
\description{
The function selects a provided list individual.
}
\usage{
select_haplotype(pop,haplotypes_ids)
}
\arguments{
\item{pop}{Required: an  input haplotype object}
\item{haplotypes_ids}{Required: a list or indices of individuals to be selected}}
\value{
A haplotype object of the selected individuals
}
\examples{
haplotype_sub1 = select_haplotype(haplotype1,c("pop1_107","pop1_167","pop1_252","pop1_288","pop1_311"))}
