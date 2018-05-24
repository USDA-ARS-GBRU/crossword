% File library/man/add_button.Rd
\name{add_button}
\alias{add_button}
\title{add_button}
\description{
adding button trigering dialog box of inputs
}
\usage{
add_button(name_function,name_output,a1,a2,e,execute_now,com)
}
\arguments{
\item{name_function}{Required: the function name}
\item{name_output}{Required: initial name of output, "" is used for function without outputs}
\item{a1}{Required: array of names of input arguments}
\item{a2}{Required: array of initial values of input arguments}
\item{e}{Required: gedit object to receive the output}
\item{execute_now}{Optional: if TRUE, the command will be execute directory instead of outputing to gedit object, default is FALSE}
\item{com}{Optional: a string holds the comments}}
\examples{
add_button("draw_haplotypes","",c("haplotypes"),c("@@pop2"),t0)}
