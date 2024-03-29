#' Collapse a matrix population model to a smaller number of stages
#'
#' Collapse a matrix population model to a smaller number of stages. For
#' instance, to compare properties of multiple projection matrices with
#' different numbers of stages, one might first collapse those matrices to a
#' standardized set of stages (e.g., propagule, pre-reproductive, reproductive,
#' and post-reproductive). The transition rates in the collapsed matrix are a
#' weighted average of the transition rates from the relevant stages of the
#' original matrix, weighted by the relative proportion of each stage class
#' expected at the stable distribution.
#'
#' @param matU The survival component of a matrix population model (i.e., a
#'   square projection matrix reflecting survival-related transitions; e.g.,
#'   progression, stasis, and retrogression)
#' @param matF The sexual component of a matrix population model (i.e., a square
#'   projection matrix reflecting transitions due to sexual reproduction)
#' @param matC The clonal component of a matrix population model (i.e., a square
#'   projection matrix reflecting transitions due to clonal reproduction).
#'   Defaults to \code{NULL}, indicating no clonal reproduction (i.e.,
#'   \code{matC} is a matrix of zeros).
#' @param collapse A list giving the mapping between stages of the original
#'   matrix and the desired stages of the collapsed matrix (e.g., \code{list(1,
#'   2:3, 4)}). Original stages may be passed as either indices or stage names
#'   corresponding to stage index or name in \code{matU}, \code{matF} and
#'   \code{matC}). Names given to the elements of \code{collapse} are used as
#'   stage names in the new, collapsed matrix.
#'
#'   See \emph{Missing Stages} for handling of \code{NA} within \code{collapse}.
#' @return A list with four elements:
#'   \item{matA}{Collapsed projection matrix}
#'   \item{matU}{Survival component of the collapsed projection matrix}
#'   \item{matF}{Sexual reproduction component of the collapsed projection
#'   matrix}
#'   \item{matC}{Clonal reproduction component of the collapsed projection
#'   matrix}
#'
#' @section Missing Stages:
#' The collapsed matrix will always be of dimension \code{length(collapse)},
#' even if one or more elements of the \code{collapse} argument is \code{NA}
#' (corresponding to a desired stage of the collapsed matrix that is not present
#' in the original matrix). In the collapsed matrix, any row/column
#' corresponding to a missing stage will be coerced to \code{NA}.
#'
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author William K. Petry <wpetry@@ncsu.edu>
#'
#' @family transformation
#'
#' @references Salguero-Gomez, R. & Plotkin, J. B. 2010. Matrix dimensions bias
#'   demographic inferences: implications for comparative plant demography. The
#'   American Naturalist 176, 710-722. <doi:10.1086/657044>
#'
#' @note This method of collapsing a matrix population model preserves the
#'   equilibrium population growth rate (\eqn{lambda}) and relative stable
#'   distribution, but is not expected to preserve other traits such as relative
#'   reproductive values, sensitivities, net reproductive rates, life
#'   expectancy, etc.
#'
#' @seealso \code{\link{mpm_standardize}}
#' @examples
#' data(mpm1)
#'
#' # check which stages reproductive
#' repro_stages(matR = mpm1$matF)
#'
#' # collapse reproductive stages (3 and 4) into single stage
#' mpm_collapse(
#'   matU = mpm1$matU, matF = mpm1$matF,
#'   collapse = list(1, 2, 3:4, 5)
#' )
#'
#' # use stage names instead, and name stages in the collapsed matrix
#' mpm_collapse(
#'   matU = mpm1$matU, matF = mpm1$matF,
#'   collapse = list(
#'     seed = "seed", vegetative = "small",
#'     flowering = c("medium", "large"),
#'     dormant = "dormant"
#'   )
#' )
#'
#' @export mpm_collapse
mpm_collapse <- function(matU, matF, matC = NULL, collapse) {
  # validate arguments
  checkValidMat(matU)
  checkValidMat(matF)
  if (!is.null(matC)) {
    checkValidMat(matC, warn_all_zero = FALSE)
  }
  checkValidStages(matU, stages = collapse)

  # populate matC with zeros, if NULL
  if (is.null(matC)) {
    outC <- FALSE
    matC <- matrix(0, nrow = nrow(matU), ncol = ncol(matU))
    if (!is.null(dimnames(matU))) {
      dimnames(matC) <- dimnames(matU)
    }
  }

  # sum components to matA
  matA <- matU + matF + matC

  # dimensions of original and collapse matrices
  originalDim <- nrow(matA)
  collapseDim <- length(collapse)

  P <- matrix(0, nrow = collapseDim, ncol = originalDim)
  if (!is.null(names(collapse))) {
    rownames(P) <- names(collapse)
  }

  # convert `collapse` names to corresponding row/col numbers if needed
  if (all(vapply(collapse, is.character, logical(1)))) {
    collapse <- lapply(collapse, function(x) which(colnames(matU) %in% x))
  }

  for (i in 1:collapseDim) {
    columns <- as.numeric(collapse[[i]])
    if (!is.na(columns[1])) {
      P[i, columns] <- 1
    }
  }

  Q <- t(P)
  w <- stable.stage(matA)

  columns <- which(colSums(Q) > 1)

  for (j in columns) {
    rows <- which(Q[, j] == 1)
    for (i in rows) {
      Q[i, j] <- w[i] / sum(w[rows])
    }
  }

  # replace missing rows/cols with NA
  if (anyNA(collapse)) {
    i <- which(is.na(collapse))
    P[i, ] <- rep(NA_real_, originalDim)
    Q[, i] <- rep(NA_real_, originalDim)
  }

  # collapse
  collapseA <- P %*% matA %*% Q
  collapseU <- P %*% matU %*% Q
  collapseF <- P %*% matF %*% Q
  collapseC <- P %*% matC %*% Q

  return(list(
    matA = collapseA,
    matU = collapseU,
    matF = collapseF,
    matC = collapseC
  ))
}
