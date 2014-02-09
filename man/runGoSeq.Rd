\name{runGoseq}
\alias{runGoSeq}
\title{runGoSeq}
\description{runGoSeq}
\usage{
runGoSeq(de.genes,all.genes,onto,genome,ID,rep=1000,p.adj="BH",adj.p.val=.05)
}
\arguments{
  \item{de.genes}{List of Differentially expressed IDs to test for GO term enrichment and depletion.}
  \item{all.genes}{List of IDs to use as the background.}
  \item{onto}{Which ontology to test, one of "GO:BP", "GO:MF" or "GO:CC".}
  \item{genome}{Which genome to use e.g. "hg19" - passed into getgo.}
  \item{ID}{Which type of ID e.g. "ensGene"- passed into getgo.}
  \item{p.adj}{Which method of p value adjustment for multiple testing to use - for details see \code{help("p.adjust")}.}
  \item{adj.p.val}{Which p value to use.}
\value{Select background explicitly for those genes annotated to root term of ontology.}
\author{Hugh French}

