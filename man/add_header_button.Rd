% File library/man/add_header_button.Rd
\name{add_header_button}
\alias{add_header_button}
\title{add_header_button}
\description{
The function create a dialog box of headers \cr
The header dialog box contains six parts: \cr
1. input/output section \cr
2. chromosome information \cr
3. founder genotypes \cr
4. phenotyping \cr
5. control buttons of submission or cancellation \cr
6. help section}
\usage{
add_header_button(pop,output,parental_genotypes)
}
\arguments{
\item{e}{Required: gedit object to receive the output}}
\examples{
add_header_button(gedit())}
