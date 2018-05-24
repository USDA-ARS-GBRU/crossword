% File library/man/count_chr_length.Rd
\name{count_chr_length}
\alias{count_chr_length}
\title{count_chr_length}
\description{
The function calculates the chromosmal lengths in bp
}
\usage{
count_chr_length(fa,output)
}
\arguments{
\item{fA}{Required: the genome fasta file}
\item{output}{Required: the output path/name}}
\value{
the genomic location.
}
\examples{
count_chr_length("peanut.fasta","chr_siz.txt")
}
