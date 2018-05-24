% File library/man/vcf2hapmap.Rd
\name{vcf2hapmap}
\alias{vcf2hapmap}
\title{vcf2hapmap}
\description{
converting vcf file format to standard hapmap format
}
\usage{
vcf2hapmap(input,output)
}
\arguments{
\item{input}{Required: the input vcf file}
\item{output}{Optional: output file name, if absent the output file will have the same name of input with an extension of .hapmap}}
\examples{
vcf2hapmap("peanut.vcf")}
