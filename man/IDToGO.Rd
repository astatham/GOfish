\name{IDToGO}
\alias{IDToGO}
\title{IDToGO}
\description{IDToGO}
\usage{
IDToGO(ids, universe=NULL, ontology=c("BP", "MF", "CC"), ID=c("entrez", "genbank", "alias", "ensembl", "symbol", "genename", "unigene"), speciesdb=c("org.Mm.eg.db", "org.Hs.eg.db"), n=NULL, p.adj=NULL, p.cutoff=0.05, orderBy=c("classic", "elim", "parentchild"), tests=NULL)
}
\arguments{
  \item{ids}{List of IDs to test for GO term enrichments}
  \item{universe}{List of IDs to use as the background. If \code{NULL} then use all IDs listed in the \code{speciesdb}}
  \item{ontology}{Which ontology to test}
  \item{ID}{Which type of ID - passed into \code{annFUN}}
  \item{speciesdb}{Which species database to use}
  \item{n}{How many GO terms to return after ordering. If 0 return all hits, if \code{NULL} n is chosen by the p value cutoff specified by \code{p.cutoff}}
  \item{p.adj}{Which method of p value adjustment for multiple testing to use - for details see \code{help("p.adjust")}. Ignored if \code{n} is specified.}
  \item{p.cutoff}{What p value cutoff to use when reporting significant GO terms. Ignored if \code{n} is specified.}
  \item{orderBy}{Which of the 3 testing algorithms to order the results in the table by. This algorithms p value is the one used as a p value cutoff when \code{n} is \code{NULL}.}
  \item{tests}{Which of the 3 testing algorithms to use (and return the results of) in addition to the one chosen in \code{orderBy}.}
}
\value{Good question}
\author{Aaron Statham}

