% File library/man/physical2genomic.Rd
\name{physical2genomic}
\alias{physical2genomic}
\title{physical2genomic}
\description{
physical2genomic function creates data structure of chosomsomes, within it, each chrosome has four data sets: 'data' contains \cr
1. "data": a data table of window information, gene frequencies and genomic loci \cr
2. "smooth": a vector of smooth function \cr
3. "size": chromsome size \cr
4. "length": chromosme length \cr
}
\usage{
physical2genomic(input,chr_stat,chr_length,window_size)
}
\arguments{
\item{input}{Required: the input gff file.}
\item{chr_stat}{Optional: an optional string argument containing a chromosome stat file which is a two-cloumn tab delimited file, \cr
the first columns is chromosome ID and the second is the size in bp.\cr
if the file was not specified, the chromsome size will be estimate based on the last gene locus by default}
\item{chr_length}{Optional: an argument containing the chromsome length information, the default is 100 CM \cr
if the argument is just number: the number will be used as a fixed value accross all chromosomes \cr
if the agrument is string: the string should has a a two-column tab delmited file, cr
the first columns is chromosome ID and the second is the length in CM.}
\item{window_size}{Optional: an optional integer for windows size, the default is 100000 bp .}}
\value{
A data structure of four data sets: "data", "smooth", "size" and "length"
}
\examples{
a1 = physical2genomic("peanut.gff")
a2 = physical2genomic("peanut.gff","chr_siz.txt")
a3 = physical2genomic("peanut.gff","chr_siz.txt",150,200000)
a4 = physical2genomic("peanut.gff","chr_siz.txt","chr_len.txt",200000)}