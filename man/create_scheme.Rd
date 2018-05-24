% File library/man/create_scheme.Rd
\name{create_scheme}
\alias{create_scheme}
\title{create_scheme}
\description{
The function creates scheme from the provided parental information: \cr
1. single parent1 X single parent2 \cr
2. single parent1 X list of parents2 \cr
3. list of parents1 X single parent2 \cr
All combinations of list of parents1 X list of parents2 is not alowed
}
\usage{
create_scheme(P1,P2,N,id)
}
\arguments{
\item{P1}{Required: parent/s 1, a name or a list of names of the parents 1}
\item{P2}{Optional: parent/s 2, a name or a list of names of the parents 2}
\item{F}{Required: number of repeated crossing generations, all crosses between founders will considered F1, the selfing will be continued for F1 up to F number of generations}
\item{N}{Required: number of individuals in each repeated generation}
\item{S}{Optional: number of lines in the segregating families descending form the last generation, F \cr
If 0, no segregant families will be produced}
\item{id}{Required: ID of the current breeding procedure}}
\value{
Output is a structured object of:
1. scheme: contains the pedigree information of all haplotypes in a data frame of four columns: \cr
	1.haplotype IDs \cr
	2.generation \cr
	3.parent 1 \cr
	4.parent 2 \cr
}
\examples{
P1 = c("Gregory","TifGurd","Tifrunner")
scheme = create_scheme(P1,"GTC20",200,'pop1')
}
