\name{ssa}
\alias{ssa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit SSA model to a time series%%  ~~function to do ... ~~
}
\description{
Fit a Singular Spectrum Analysis (SSA) model to a univariate time series. This function provides options for interactively choosing the number of components to retain, controlling verbosity of the output, and producing diagnostic and other plots.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ssa(tsdata, window, train.periods,
    q = NULL, prop.var.exp = NULL,
    verbose = FALSE, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tsdata}{
  A vector containing univariate time series data.
%%     ~~Describe \code{tsdata} here~~
}
  \item{window}{
  Window length for the SSA model. Should be a value between 2 and half the length of the time series.
%%     ~~Describe \code{window} here~~
}
  \item{train.periods}{
  Number of data points (at the start of the time series) used for training the SSA model.
%%     ~~Describe \code{train.periods} here~~
}
  \item{q}{
  Number of SSA components to retain for modelling purposes, if known in advance. Ignored if \code{prop.var.exp} is specified.
%%     ~~Describe \code{q} here~~
}
  \item{prop.var.exp}{
  Required minimum proportion of variation accounted for by the SSA model. Should be a value in the range (0, 1].
%%     ~~Describe \code{prop.var.exp} here~~
}
  \item{verbose}{
  Logical (\code{TRUE}/\code{FALSE}). Should modelling statistics be shown in the console?
%%     ~~Describe \code{verbose} here~~
}
  \item{plot}{
  Logical (\code{TRUE}/\code{FALSE}). Should modelling related plots be shown?
%%     ~~Describe \code{plot} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list with the components:
\item{tsdata}{Original time series data.}
\item{fitted}{Fitted/predicted values from the SSA model.}
\item{values}{Singular values.}
\item{vectors}{Singular vectors.}
\item{q}{Number of components retained in SSA model.}
\item{projmat}{Projection matrix for the SSA model.}
\item{rowmeans}{Row means of the training data trajectory matrix.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Theo Pepler
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x <- rpois(n = 120, lambda = 20)
ssa(tsdata = x, window = 12, train.periods = 60, verbose = TRUE, plot = TRUE)

# NOT RUN --- Example where proportion variation is specified
# ssa(tsdata = x, window = 12, train.periods = 60, prop.var.exp = 0.8, verbose = TRUE, plot = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
