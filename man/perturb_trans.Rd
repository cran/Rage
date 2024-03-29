% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perturb_trans.R
\name{perturb_trans}
\alias{perturb_trans}
\title{Perturbation analysis of transition types within a matrix population model}
\usage{
perturb_trans(
  matU,
  matF,
  matC = NULL,
  posU = matU > 0,
  posF = matF > 0,
  posC = matC > 0,
  exclude_row = NULL,
  exclude_col = NULL,
  pert = 1e-06,
  type = "sensitivity",
  demog_stat = "lambda",
  ...
)
}
\arguments{
\item{matU}{The survival component submatrix of a MPM (i.e., a square
projection matrix reflecting survival-related transitions; e.g.,
progression, stasis, and retrogression).}

\item{matF}{The sexual component submatrix of a MPM (i.e., a square
projection matrix reflecting transitions due to sexual reproduction).}

\item{matC}{The clonal component submatrix of a MPM (i.e., a square
projection matrix reflecting transitions due to clonal reproduction).
Defaults to \code{NULL}, indicating no clonal reproduction possible.}

\item{posU}{A logical matrix of the same dimension as \code{matU}, with
elements indicating whether a given \code{matU} transition is possible
(\code{TRUE}) or not (\code{FALSE}). Defaults to \code{matU > 0} (see
Details).}

\item{posF}{A logical matrix of the same dimension as \code{matF}, with
elements indicating whether a given \code{matF} transition is possible
(\code{TRUE}) or not (\code{FALSE}). Defaults to \code{matF > 0} (see
Details).}

\item{posC}{A logical matrix of the same dimension as \code{matC}, with
elements indicating whether a given \code{matC} transition is possible
(\code{TRUE}) or not (\code{FALSE}). Defaults to \code{matC > 0} (see
Details).}

\item{exclude_row}{A vector of row indices or stages names indicating stages
for which transitions \emph{to} the stage should be excluded from
perturbation analysis. Alternatively, a logical vector of length
\code{nrow(matU)} indicating which stages to include \code{TRUE} or exclude
\code{FALSE} from the calculation. See section \emph{Excluding stages}.}

\item{exclude_col}{A vector of column indices or stages names indicating
stages for which transitions \emph{to} the stage should be excluded from
perturbation analysis. Alternatively, a logical vector of length
\code{ncol(matU)} indicating which stages to include \code{TRUE} or exclude
\code{FALSE} from the calculation. See section \emph{Excluding stages}.}

\item{pert}{The magnitude of the perturbation (defaults to \code{1e-6}).}

\item{type}{An argument defining whether to return `sensitivity` or
`elasticity` values. Defaults to `sensitivity`.}

\item{demog_stat}{An argument defining which demographic statistic should be
used, as in "the sensitivity/elasticity of \code{demog_stat} to matrix
element perturbations." Defaults to the per-capita population growth rate
at equilibrium (\eqn{lambda}). Also accepts a user-supplied function that
performs a calculation on a MPM and returns a single numeric value.}

\item{...}{Additional arguments passed to the function \code{demog_stat}.}
}
\value{
A list with 5 elements: \item{stasis}{The sensitivity or elasticity
  of \code{demog_stat} to stasis.} \item{retrogression}{The sensitivity or
  elasticity of \code{demog_stat} to retrogression.} \item{progression}{The
  sensitivity or elasticity of \code{demog_stat} to progression.}
  \item{fecundity}{The sensitivity or elasticity of \code{demog_stat} to
  sexual fecundity.} \item{clonality}{The sensitivity or elasticity of
  \code{demog_stat} to clonality.}
}
\description{
Calculates the summed sensitivities or elasticities for various
transition types within a matrix population model (MPM), including stasis,
retrogression, progression, fecundity, and clonality.

Sensitivities or elasticities are calculated by perturbing elements of the
MPM and measuring the response of the per-capita population growth rate at
equilibrium (\eqn{\lambda}), or, with a user-supplied function, any other
demographic statistic.
}
\details{
A transition rate of \code{0} within a matrix population model can
either indicate that the transition is not possible in the given life cycle
(e.g., tadpoles never revert to eggs), or that the transition is possible but
was estimated to be \code{0} in the relevant population and time period.
Because transition rates of zero \emph{do} generally yield non-zero
sensitivities, it is important to distinguish between structural (i.e.
impossible) zeros and sampled zeros when summing multiple sensitivities for a
given process (e.g., progression/growth).

By default, the \code{perturb_} functions assume that a transition rate of
\code{0} indicates an impossible transition, in which case the sensitivity
for that transition will not be included in any calculation. Specifically,
the arguments \code{posX} are specified by the logical expression \code{(matX
> 0)}. If the matrix population model includes transitions that are possible
but estimated to be \code{0}, users should specify the \code{posX}
argument(s) manually.

If there are no possible transitions for a given process (e.g., clonality, in
many species), the value of sensitivity or elasticity returned for that
process will be \code{NA}.
}
\section{Excluding stages}{
 It may be desirable to exclude one or more stages
  from the calculation. For instance, we might not believe that 'progression'
  to a dormant stage class truly reflects progression. In this case we could
  exclude transitions \emph{to} the dormant stage class using the argument
  \code{exclude_row}. We may or may not want to ignore progression
  transitions \emph{from} the dormant stage class, which can be done in a
  similar way using the argument \code{exclude_col}. The \code{exclude_}
  arguments simply set the relevant row or column of the \code{posX}
  arguments to \code{FALSE}, to prevent those transitions from being used in
  subsequent calculations.
}

\examples{
matU <- rbind(
  c(0.1, 0, 0, 0),
  c(0.5, 0.2, 0.1, 0),
  c(0, 0.3, 0.3, 0.1),
  c(0, 0, 0.5, 0.6)
)

matF <- rbind(
  c(0, 0, 1.1, 1.6),
  c(0, 0, 0.8, 0.4),
  c(0, 0, 0, 0),
  c(0, 0, 0, 0)
)


perturb_trans(matU, matF)

# Use a larger perturbation than the default of 1e-6.
perturb_trans(matU, matF, pert = 0.01)

# Calculate the sensitivity/elasticity of the damping ratio to perturbations.
# First, define function for damping ratio:
damping <- function(matA) {
  eig <- eigen(matA)$values
  dm <- rle(Mod(eig))$values
  return(dm[1] / dm[2])
}

# Second, run the perturbation analysis using demog_stat = "damping".
perturb_trans(matU, matF, demog_stat = "damping")

}
\seealso{
Other perturbation analysis: 
\code{\link{perturb_matrix}()},
\code{\link{perturb_stochastic}()},
\code{\link{perturb_vr}()},
\code{\link{pop_vectors}()}
}
\author{
Rob Salguero-Gómez <rob.salguero@zoo.ox.ac.uk>

Patrick Barks <patrick.barks@gmail.com>
}
\concept{perturbation analysis}
