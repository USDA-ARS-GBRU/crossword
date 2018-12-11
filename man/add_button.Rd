% File library/man/add_button.Rd
\name{add_button}
\alias{add_button}
\title{add_button}
\description{
Adding button that trigger a customized dialog box. \cr
The dialog box contains three parts: \cr
1. An optional single output object space \cr
2. Single or multiple input objects \cr
3. control buttons of submission or cancellation}
\usage{
add_button(name_function,name_output,a1,a2,e,execute_now,com,typein,buttons)
}
\arguments{
\item{name_function}{Required: the function name}
\item{name_output}{Required: initial name of output, "" is used for function without outputs}
\item{a1}{Required: array of names of input arguments}
\item{a2}{Required: array of initial values of input arguments}
\item{e}{Required: gedit object to receive the output}
\item{execute_now}{Optional: if TRUE, the command will be execute directory instead of outputting to gedit object, default is FALSE}
\item{com}{Optional: a string holds the comments}
\item{typein}{Optional: types of the entry box, defualt is "gedit" object}
\item{buttons}{Optional: string object is used if some of the entries require extra dialog box \cr
the user may select from "tkgetSaveFile", "tkgetOpenFile", "tkchooseDirectory", OR NA}}
\examples{
add_button("draw_haplotypes","",c("haplotypes"),c("@@pop2"),t0)}
