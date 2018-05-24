% File library/man/draw_haplotypes.Rd
\name{draw_haplotypes}
\alias{draw_haplotypes}
\title{draw_haplotypes}
\description{
The function plots haplotypes
}
\usage{
draw_haplotypes(haplotypes,output_folder,parental_genotypes,heterozygous,by_chromosomes,im_type)
}
\arguments{
\item{haplotypes}{Required: a haplotype object input.}
\item{output_folder}{Required: an output folder to be created and store the images.}
\item{parental_genotypes}{Required: the parental genotyping object}
\item{heterozygous}{Optional: a logical argument to provide if the crop is heterozygous, default: TRUE}
\item{by_chromosomes}{Optional: a logical argument to provide if the haplotypes will be sorted in the graphs by chromsomes, default: TRUE}
\item{im_type}{Optional: the type of images}}

\examples{
founders = read.table("parents.txt")[,1]
draw_haplotypes(haplotypes,'output_haplo',parental_genotypes)}
