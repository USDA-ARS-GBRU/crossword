% File library/man/create_genome.Rd
\name{create_genome}
\alias{create_genome}
\title{create_genome}
\description{
create_genome creates two sets of genomes for every individual \cr
It takes a fasta genome file and genotype object as hapmap file. \cr
For every individual in the genotype file, two genome fasta files will be produced for the two possible alleles in for every loci.}
\usage{
create_genome(fa,genotype,output,haploid_only,max_total_size_in_gb)
}
\arguments{
\item{fa}{Required: the genome fasta file}
\item{genotype}{Required: a hapmap file}
\item{output}{Required: a prefix of output files}
\item{haploid_only}{Optional: if haploid or diploid genome will be produced, default is TRUE}
\item{max_total_size_in_gb}{Optional: size quota which should not be exceeded, default is 10 gb}}
\value{
the genomic location.
}
\examples{
fa = "Gmax_275_v2.0.fa"
genotype = "pop8_marker_set.hapmap.hapmap"
output="genome"
create_genome(fa,genotype,output)}
