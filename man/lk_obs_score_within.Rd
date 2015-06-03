\name{lk_obs_score_within}
\alias{lk_obs_score_within}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute observed log-likelihood and score}
\description{Function used within \code{est_multi_poly} to compute observed log-likelihood and score.}
\usage{lk_obs_score_within(par_comp,lde1,lde2,lpar,lga1,lga2,S,R,yv,k1,k2,rm1,rm2,lv,J,
	                       fv,link,disc,indga1,indga2,glob,refitem,miss,ltype,XX1dis,
	                       XX2dis,Xlabel,ZZ0,fort)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par_comp}{complete vector of parameters}
  \item{lde1}{length of de1}
  \item{lde2}{length of de2}
  \item{lpar}{length of par}
  \item{lga1}{length of ga1}
  \item{lga2}{length of ga2}
  \item{S}{matrix of responses}
  \item{R}{matrix of observed responses indicator}
  \item{yv}{vector of frequencies}
  \item{k1}{number of latent classes for the 1st latent variable}
  \item{k2}{number of latent classes for the 2nd latent variable}
  \item{rm1}{number of dimensions for the 1st latent variable}
  \item{rm2}{number of dimensions for the 2nd latent variable}
  \item{lv}{number of response categories for each item}
  \item{J}{number of items}
  \item{fv}{indicator of constrained parameters}
  \item{link}{link function}
  \item{disc}{presence of discrimination parameter}
  \item{indga1}{indicator of ga1 parameters}
  \item{indga2}{indicator of ga2 parameters}
  \item{glob}{indicator of global parametrization for the covariates}
  \item{refitem}{vector of reference items}
  \item{miss}{indicator of presence of missing responses}
  \item{ltype}{type of logit}
  \item{XX1dis}{array of covariates for the 1st latent variable}
  \item{XX2dis}{array of covariates for the 2nd latent variable}
  \item{Xlabel}{indicator for covariate configuration}
  \item{ZZ0}{design matrix}
  \item{fort}{to use Fortran}
}

%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}

\value{
\item{lk}{log-likelihood function}
\item{sc}{score vector}
}

\author{Francesco Bartolucci - University of Perugia (IT)}

%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
