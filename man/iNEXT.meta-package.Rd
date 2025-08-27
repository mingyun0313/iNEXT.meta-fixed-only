\name{iNEXT.meta-package}
\alias{iNEXT.meta-package}
\alias{iNEXT.meta}
\docType{package}
\title{
  Meta-analysis of the difference between two treatments in interpolation and extrapolation for beta diversity across three dimensions
}
\description{
  \code{iNEXT.meta} (Meta-analysis of the difference between two treatments in interpolation and extrapolation for beta diversity across three dimensions) is an R package that extends the concepts of \code{iNEXT.beta3D} (Chao et al., 2023) to meta-analysis (fixed-effect model). \cr

  \code{iNEXT.meta} is primarily designed to integrate results comparing diversity differences between two treatments across multiple studies. It enables the calculation of three dimensions (3D) of diversity: taxonomic diversity (TD), phylogenetic diversity (PD) and functional diversity (FD), at a specified level of sample coverage. \cr

  This package provides three main functions: \cr

  1. \code{DataInfobeta3Dmeta} provides basic data information in each combination of site and treatment for (1) the gamma reference sample in the pooled assemblage, and (2) the alpha reference sample in the joint assemblage. \cr

  2. \code{iNEXTbeta3Dmeta} estimates the difference of standardized 3D diversity with common sample coverage (for alpha, beta, gamma diversity, and four classes of dissimilarity measures) between two treatments, and fit a fixed-effect model to perform meta-analysis. \cr

  3. \code{ggiNEXTmeta} visualizes the output from the function \code{iNEXTbeta3Dmeta} by providing forest plot of the difference between two treatments of standardized 3D diversity in each study/site and meta-analysis. \cr

}

\keyword{ package }
