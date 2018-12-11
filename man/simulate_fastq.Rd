% File library/man/simulate_fastq.Rd
\name{simulate_fastq}
\alias{simulate_fastq}
\title{simulate_fastq}
\description{
The function simulates illumina fastq files (a paired-end or single-end) using ART.}
\usage{
simulate_fastq(single_reads,input,read_len,fold,art_binary_location,genotypes)
}
\arguments{
\item{single_reads}{Optional: it direct the function to use paired-end, default FALSE, or single-end, TRUE, reads simulation}
\item{genotypes}{Optional: it should be a hapmap file. If supplied, the fastq files will be simulated from real_time genomes produced from input fasta and genotypes containing in this hapmap file}
\item{input}{Required: a fasta file or a directory contains a group of fasta files to be simulated, in this case the fold will be devided by number of the files in such directory}
\item{read_len}{Optional: estimated read length, default is 250}
\item{fold}{Optional: estimated read length, default is 30x}
\item{art_binary_location}{Optional: the location of ART binary softwar, , default is art version implemented with this R package, 2016.06.05}
}
\examples{
simulate_fastq(FALSE,"test/a1.fa",read_len=250,fold=1.32"~/Documents/apps/art_bin_GreatSmokyMountains")
simulate_fastq(FALSE,"test",read_len=250,fold=1.32,"~/Documents/apps/art_bin_GreatSmokyMountains")}
