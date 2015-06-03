\name{est_multi_poly_within}
\alias{est_multi_poly_within}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Estimate LC IRT model for dichotomous and polytomous responses under within multidimensionality}

\description{The function performs maximum likelihood estimation of the parameters of the IRT models
assuming a discrete distribution for the ability and within multidimensionality. Every ability level corresponds to a latent
class of subjects in the reference population. Maximum likelihood estimation is based on Expectation-
Maximization algorithm.}

\usage{est_multi_poly_within(S, yv = rep(1,ns), k1, k2, X = NULL, start = 0, link = 1,
        disc = 0, difl = 0, multi1 = 1:J, multi2 = 1:J, piv1 = NULL, piv2 = NULL, 
        Phi = NULL, gac = NULL, De = NULL, fort = FALSE, tol = 10^-10, disp = FALSE,
        output = FALSE, out_se = FALSE, glob = FALSE)}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{S}{matrix of all response sequences observed at least once in the sample and listed row-by-row
           (use NA for missing response)}
  \item{yv}{vector of the frequencies of every response configuration in \code{S}}
  \item{k1}{number of ability levels (or latent classes) for the 1st latent variable}
  \item{k2}{number of ability levels (or latent classes) for the 2nd latent variable}
  \item{X}{matrix of covariates that affects the weights}
  \item{start}{method of initialization of the algorithm (0 = deterministic, 1 = random, 2 = arguments given as
               input)}
  \item{link}{type of link function (0 = no link function, 1 = global logits, 2 = local logits);
              with no link function the Latent Class model results; with global logits the Graded Response model
              results; with local logits the  Partial Credit results (with dichotomous responses, global logits
              is the same as using local logits resulting in the Rasch or the 2PL model depending on the value
              assigned to disc)}
  \item{disc}{indicator of constraints on the discriminating indices (0 = all equal to one, 1 = free)}
  \item{difl}{indicator of constraints on the difficulty levels (0 = free, 1 = rating scale parameterization)}
  \item{multi1}{matrix with a number of rows equal to the number of dimensions and elements in each row
               equal to the indices of the items measuring the dimension corresponding to that row for the 1st latent variable}
  \item{multi2}{matrix with a number of rows equal to the number of dimensions and elements in each row
               equal to the indices of the items measuring the dimension corresponding to that row for the 2nd latent variable}
  \item{piv1}{initial value of the vector of weights of the latent classes (if start=2) for the 1st latent variable}
  \item{piv2}{initial value of the vector of weights of the latent classes (if start=2) for the 2nd latent variable}
  \item{Phi}{initial value of the matrix of the conditional response probabilities (if start=2)}
  \item{gac}{initial value of the complete vector of discriminating indices (if start=2)}
  \item{De}{initial value of regression coefficients for the covariates (if start=2)}
  \item{fort}{to use Fortran routines when possible}
  \item{tol}{tolerance level for checking convergence of the algorithm as relative difference between 
  	         consecutive log-likelihoods}
  \item{disp}{to display the likelihood evolution step by step}
  \item{output}{to return additional outputs (Phi,Pp,Piv)}
  \item{out_se}{to return standard errors}
  \item{glob}{to use global logits in the covariates}  
}

%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}

\value{
\item{piv1}{estimated vector of weights of the latent classes (average of the weights in case of model with covariates) for the 1st latent variable}
\item{piv2}{estimated vector of weights of the latent classes (average of the weights in case of model with covariates) for the 2nd latent variable}
\item{Th1}{estimated matrix of ability levels for each dimension and latent class for the 2nd latent variable}
\item{Th2}{estimated matrix of ability levels for each dimension and latent class for the 2nd latent variable}
\item{Bec}{estimated vector of difficulty levels for every item (split in two vectors if difl=1)}
\item{ga1c}{estimated vector of discriminating indices for every item (with all elements equal to 1
	       with Rasch parametrization) for the 1st latent variable}
\item{ga2c}{estimated vector of discriminating indices for every item (with all elements equal to 1
	       with Rasch parametrization) for the 2nd latent variable}
\item{fv1}{vector indicating the reference item chosen for each latent dimension for the 1st latent variable}
\item{fv2}{vector indicating the reference item chosen for each latent dimension for the 2nd latent variable}
\item{Phi}{array of the conditional response probabilities for every item and each of the k1*k2 latent classes}
\item{De}{matrix of regression coefficients for the multinomial logit model on the class weights}
\item{Piv}{matrix of the weights for every response configuration (if output=TRUE)}
\item{Pp}{matrix of the posterior probabilities for each response configuration and latent class
	      (if output=TRUE)}
\item{lk}{log-likelhood at convergence of the EM algorithm}
\item{np}{number of free parameters}
\item{aic}{Akaike Information Criterion index}
\item{bic}{Bayesian Information Criterion index}
\item{ent}{Etropy index to measure the separation of classes}
\item{lkv}{Vector to trace the log-likelihood evolution across iterations (if output=TRUE)}
\item{seDe1}{Standard errors for De (if output=TRUE) for the 1st latent variable}
\item{seDe2}{Standard errors for De (if output=TRUE) for the 2nd latent variable}
\item{separ}{Standard errors for vector of parameters containing Th and Be (if out_se=TRUE)}
\item{sega1}{Standard errors for vector of discrimination indices (if out_se=TRUE) for the 1st latent variable}
\item{sega2}{Standard errors for vector of discrimination indices (if out_se=TRUE) for the 2nd latent variable}
\item{Vn}{Estimated variance-covariance matrix for all parameter estimates (if output=TRUE)}
}

\references{
Bartolucci, F. (2007), A class of multidimensional IRT models for testing unidimensionality and clustering
items, \emph{Psychometrika}, \bold{72}, 141-157.

Bacci, S., Bartolucci, F. and Gnaldi, M. (2014), A class of Multidimensional Latent Class IRT models for 
ordinal polytomous item responses, \emph{Communication in Statistics - Theory and Methods}, \bold{43}, 787-800.
}

\author{Francesco Bartolucci, Silvia Bacci - University of Perugia (IT)}

%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{

\dontrun{
## Fit the model under different situations for NAEP data
data(naep)
multi1 = 1:12
multi2 = 2:12
# Rasch model and 2PL with two dimensions (within)
out10 = est_multi_poly_within(naep,k1=3,k2=3,multi1=multi1,multi2=multi2,link=1,disp=TRUE,
                  out_se=TRUE) 
out11 = est_multi_poly_within(naep,k1=3,k2=3,multi1=multi1,multi2=multi2,link=1,disc=1,
                  disp=TRUE,out_se=TRUE) 
}

\dontrun{
# fit the model under different situations for RLMS data
# example of use of the function to account for non-ignorable missing responses
data(RLMS)
X = RLMS[,1:4]
Y = RLMS[,6:9]
YR = cbind(Y,1*(!is.na(Y)))
multi1 = 1:4
multi2 = 5:8
# MAR model
out0 = est_multi_poly_within(YR,X=X,k1=3,k2=2,multi1=multi1,multi2=multi2,link=1,difl=1,
                 disc=1,glob=1,disp=TRUE,out_se=TRUE) 
# NMAR model
multi1 = 1:8
out1 = est_multi_poly_within(YR,X=X,k1=3,k2=2,multi1=multi1,multi2=multi2,link=1,difl=1,
                 disc=1,glob=1,disp=TRUE,out_se=TRUE) 
# testing effect of the latent trait on missingness
c(out0$bic,out1$bic)
test1 = out1$ga1c[-1]/out1$sega1
}
}

%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
\keyword{maximum likelihood estimation}
\keyword{Expectation-Maximization algorithm}% __ONLY ONE__ keyword per line
