% File library/man/filter_parents.Rd
\name{filter_parents}
\alias{filter_parents}
\title{filter_parents}
\description{
filter vcf or hapmap files to ceratin parents
}
\usage{
filter_parents(input,output,selected_parents_file)
}
\arguments{
\item{input}{Required: an input vcf of hapmap files}
\item{output}{Required: output file path/name}
\item{selected_parents_file}{Required: an input file of a list of parents}}
\examples{
filter_parents("parental_genotypes.hapmap","filtered_parental_genotypes.hapmap","selected_parents_file")}
