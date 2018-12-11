% File library/man/get_level.Rd
\name{get_level}
\alias{get_level}
\title{get_level}
\description{
The function returns the level factor on the input array and using the supplied level.}
\usage{
get_level(m1,level,return_list,print_warnings)
}
\arguments{
\item{m1}{Required: a scheme matrix}
\item{level}{Required: the requested level: individual, family, cross or population}
\item{return_list}{Optional: a logical argument to show if the return will be one dimension array or a list}
\item{print_warnings}{Optional: TRUE to print the warnings}}
\value{
a two-column tab delimited phenotyping object; id and value
}
\examples{
m2 = get_level(m1,"cross",FALSE)}
