% File library/man/advance.Rd
\name{advance}
\alias{advance}
\title{advance}
\description{
The function create advanced generations from a cross population
}
\usage{
advance(pop,F,outcross,level,clevel,id,CHR)
}
\arguments{
\item{pop}{Required: a cross population}
\item{F}{Required: Number of F generations}
\item{id}{Required: the new population ID}
\item{outcross}{Required: outcross level within the crop}
\item{level}{Required: the level of cross}
\item{clevel}{Required: the level of outcrosses}
\item{CHR}{Required: a two-column data frame of chromosomes' IDs and length}
}
\value{
A haplotype object (population)
}
\examples{
pop2 = advance(pop=pop1,F=5,outcross=0.1,level=0,clevel=0,id="pop2",CHR="chr_len.txt")}

