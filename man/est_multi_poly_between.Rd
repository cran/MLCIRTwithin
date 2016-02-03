\name{est_multi_poly_between}
\alias{est_multi_poly_between}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Estimate latent class item response theory (LC-IRT)  models for dichotomous and polytomous responses under between-item multidimensionality}

\description{The function performs maximum likelihood estimation of the parameters of the IRT models
assuming a discrete distribution for the ability and between-item multidimensionality. Every ability level corresponds to a latent
class of subjects in the reference population. Maximum likelihood estimation is based on Expectation-
Maximization algorithm.}

\usage{est_multi_poly_between(S, yv = rep(1,ns), k, X = NULL,
                       start = c("deterministic","random","external"),
                       link = c("global","local"), disc = FALSE,
                       difl = FALSE, multi = 1:J, piv = NULL,
                       Phi = NULL, gac = NULL, De = NULL, fort = FALSE,
                       tol = 10^-10, disp = FALSE, output = FALSE,
                       out_se = FALSE, glob = FALSE)}
        
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{S}{matrix of all response sequences observed at least once in the sample and listed row-by-row
           (use NA for missing responses)}
  \item{yv}{vector of the frequencies of every response configuration in \code{S}}
  \item{k}{number of ability levels (or latent classes) for the latent variable}
  \item{X}{matrix of covariates that affects the weights}
  \item{start}{method of initialization of the algorithm)}
  \item{link}{type of link function ("global" for global logits, "local" for local logits);
               with global logits a graded response model results; with local logits a  partial credit model results (with dichotomous responses, global  				logits is the same as using local logits resulting in the Rasch or the 2PL model depending on the value
              assigned to disc)}
  \item{disc}{indicator of constraints on the discriminating indices (FALSE = all equal to one, TRUE = free)}
  \item{difl}{indicator of constraints on the difficulty levels (FALSE = free, TRUE = rating scale parametrization)}
  \item{multi}{matrix with a number of rows equal to the number of dimensions and elements in each row
               equal to the indices of the items measuring the dimension corresponding to that row for the latent variable}
  \item{piv}{initial value of the vector of weights of the latent classes (if start="external") for the latent variable}
  \item{Phi}{initial value of the matrix of the conditional response probabilities (if start="external")}
  \item{gac}{initial value of the complete vector of discriminating indices (if start="external")}
  \item{De}{initial value of regression coefficients for the covariates (if start="external")}
  \item{fort}{to use Fortran routines when possible}
  \item{tol}{tolerance level for checking convergence of the algorithm as relative difference between 
  	         consecutive log-likelihoods}
  \item{disp}{to display the likelihood evolution step by step}
  \item{output}{to return additional outputs (Piv,Pp,lkv)}
  \item{out_se}{to return standard errors}
  \item{glob}{to use global logits in the covariates}  
}

%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}

\value{
\item{piv}{estimated vector of weights of the latent classes (average of the weights in case of model with covariates)}
\item{Th}{estimated matrix of ability levels for each dimension and latent class of the latent variable}
\item{Bec}{estimated vector of difficulty levels for every item (split in two vectors if difl=TRUE)}
\item{gac}{estimated vector of discriminating indices for every item (with all elements equal to 1
	       with Rasch parametrization)}
\item{fv}{vector indicating the reference item chosen for each latent dimension of the latent variable}
\item{Phi}{array of the conditional response probabilities for every item and each of the k latent classes}
\item{De}{matrix of regression coefficients for the multinomial logit model (or global logit model if glob=TRUE) on the class weights}
\item{Piv}{matrix of the weights for every response configuration (if output=TRUE)}
\item{Pp}{matrix of the posterior probabilities for each response configuration and latent class
	      (if output=TRUE)}
\item{lk}{log-likelhood at convergence of the EM algorithm}
\item{np}{number of free parameters}
\item{aic}{Akaike Information Criterion index}
\item{bic}{Bayesian Information Criterion index}
\item{ent}{Entropy index to measure the separation of classes}
\item{lkv}{Vector to trace the log-likelihood evolution across iterations (if output=TRUE)}
\item{seDe}{Standard errors for De (if out_se=TRUE)}  
\item{seTh}{Standard errors for vector Th (if out_se=TRUE)}
\item{seBec}{Standard errors for vector Bec of difficulty parameters (if out_se=TRUE)}
\item{sega}{Standard errors for vector gac of discrimination indices (if out_se=TRUE)}
\item{Vn}{Estimated variance-covariance matrix for all parameter estimates (if out_se=TRUE)}
}

\references{
Bartolucci, F. (2007), A class of multidimensional IRT models for testing unidimensionality and clustering
items, \emph{Psychometrika}, \bold{72}, 141-157.

Bacci, S., Bartolucci, F. and Gnaldi, M. (2014), A class of Multidimensional Latent Class IRT models for 
ordinal polytomous item responses, \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 787-800.
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
# Fit a Graded response model with two latent variables (free discrimination
# and difficulty parameters; two latent classes):
data(SF12_nomiss)
S = SF12_nomiss[,1:12]
X = SF12_nomiss[,13]
multi0 = rbind(c(1:5, 8), c(6:7,9:12))
out1 =  est_multi_poly_between(S=S,k=2,X=X,link="global",disc=TRUE,
                               multi=multi0,fort=TRUE,disp=TRUE,out_se=TRUE) 

# Display output:
out1$lk
out1$Th
out1$piv
out1$De
}

}

%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
\keyword{maximum likelihood estimation}
\keyword{Expectation-Maximization algorithm}% __ONLY ONE__ keyword per line
