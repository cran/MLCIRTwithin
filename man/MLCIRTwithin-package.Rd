\name{MLCIRTwithin-package}

\alias{MLCIRTwithin-package}

\alias{MLCIRTwithin}

\docType{package}

\title{Latent Class (LC) Item Response Theory (IRT) Models Under 'Within-Item Multi-Dimensionality'}

\description{This package provides a flexible framework for the Item Response Theory (IRT) analysis 
of dichotomous and ordinal polytomous outcomes 
under the assumption of within-item multi-dimensionality and discreteness of latent traits (abilities). 
Every level of the abilities identify a latent class of subjects. The fitting 
algorithms are based on the Expectation-Maximization (EM) paradigm and 
allow for missing responses and for different item parameterizations. The package also allows for the inclusion
individual covariates affecting the class weights.}

\details{
\tabular{ll}{
Package: \tab MultiLCIRT\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-06-02\cr
License: \tab GPL (>= 2)\cr
}
Function \code{est_multi_poly_with} performs the parameter estimation of the same model considered in the \code{R} package \code{MultiLCIRT} under withing multidimensionality.
}
\author{Francesco Bartolucci, Silvia Bacci, Michela Gnaldi - University of Perugia (IT)

Maintainer: Francesco Bartolucci <bart@stat.unipg.it>
}
% \references{
% Bartolucci, F. (2007), A class of multidimensional IRT models for testing unidimensionality and clustering
% items, \emph{Psychometrika}, \bold{72}, 141-157.

% Bacci, S., Bartolucci, F. and Gnaldi, M. (2014), A class of Multidimensional Latent Class IRT models for 
% ordinal polytomous item responses, \emph{Communication in Statistics - Theory and Methods}, \bold{43}, 787-800.

% Bartolucci, F., Bacci, S. and Gnaldi, M. (2014), MultiLCIRT: An R package for multidimensional latent class item response models, \emph{Computational Statistics and Data Analysis}, \bold{71}, 971-985.
% }
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
% \examples{
% ## Estimation of different Multidimensional LC IRT models with binary responses
% # Aggregate data
% data(naep)
% X = as.matrix(naep)
% out = aggr_data(X)
% S = out$data_dis
% yv = out$freq
% # Define matrix to allocate each item on one dimension 
% multi1 = rbind(c(1,2,9,10),c(3,5,8,11),c(4,6,7,12))
% # Three-dimensional LC Rasch model with 4 latent classes
% # less severe tolerance level to check convergence (to be modified)
% out1 = est_multi_poly(S,yv,k=4,start=0,link=1,multi=multi1,tol=10^-6)
% }
