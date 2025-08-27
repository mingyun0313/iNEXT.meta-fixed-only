\name{Bird_distM}
\alias{Bird_distM}
\docType{data}
\title{
  Species pairwise distance matrix for bird species
}
\description{
  This dataset includes a bird species pairwise symmetric distance matrix for 72 species listed in the \code{Bird_incidence_data} dataset. Each element in the matrix represents a Gower distance between two species computed from species traits. \cr
}
\usage{data("Bird_distM")}
\format{
  \code{Bird_distM} is a 72 x 72 (species by species) symmetric data.frame. Each element of the data.frame is between zero (for species pairs with identical traits) and one. \cr
}

\keyword{datasets}
