\name{create450kbg}
\alias{create450kbg}
\title{create450kbg}
\description{create450kbg}
\usage{
create450kbg(probes, groups=NULL)
}
\arguments{
  \item{probes}{Character vector of Illumina 450k IDs (e.g. "cg00000108"") constituting background}
  \item{groups}{Character vector of one or combination of "TSS1500", "TSS200", "Body", "1stExon", "5'UTR" and "3'UTR"}
}
\value{A list object of class "background" with 3 elements: "background" is a gene-symbol named vector of integers describing how many 450k probes in the constructed background are associated with that gene symbol; "probes" is a character vector of Illumina probe IDs constituting the background; "groups" is the argument passed to \code{groups}}
\author{Tim J. Peters <Tim.Peters@csiro.au>}
\examples{
\dontrun{
#For example, creating a background of promoter probes only
bg <- create450kbg(my.beta.matrix, groups=c("TSS1500", "TSS200"))}
}

