\name{regionsToIDs}
\alias{regionsToIDs}
\title{regionsToIDs}
\description{regionsToIDs}
\usage{
regionsToIDs(x, txdb, where=c("TSS", "gene"), within=100)
}
\arguments{
  \item{x}{A \code{GRanges} of regions of which to obtain the ID of genes contained within}
  \item{txdb}{A \code{TranscriptDb} made using \code{makeTranscriptDb} from the \code{GenomicFeatures} package}
  \item{where}{If "TSS", only genes whose transcription start site overlaps the regions in \code{x} are returned. If "gene" the whole gene is considered, not just the TSS.}
  \item{within}{If \code{where} is "TSS" then the number of basepairs around each TSS a region from \code{x} may overlap to be returned. Ignored if \code{where} is "gene".}
}
\value{A character vector containing gene IDs of genes contained within the supplied regions}
\author{Aaron Statham}

