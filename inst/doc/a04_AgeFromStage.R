## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(Rage) # load Rage
data(mpm1) # load data object 'mpm1'
mpm1 # display the contents

## -----------------------------------------------------------------------------
dimnames(mpm1$matU)

# extract U and F matrices
mat_U <- mpm1$matU
mat_F <- mpm1$matF

# calculate lx
lx <- mpm_to_lx(mat_U, start = 1, xmax = 30)

# calculate mx
mx <- mpm_to_mx(mat_U, mat_F, start = 1, xmax = 30)

## -----------------------------------------------------------------------------
# calculate lx
lx <- mpm_to_lx(mat_U, start = 1, xmax = 30)

# calculate mx
mx <- mpm_to_mx(mat_U, mat_F, start = 1, xmax = 30)

## ---- fig.width = 6, fig.height = 4-------------------------------------------
plot(lx, ylim = c(0, 1), type = "l", xlab = "Age")
plot(mx, type = "l", xlab = "Age")

## -----------------------------------------------------------------------------
library(Rcompadre)
data(Compadre)

# In older versions of Com(p)adre the ProjectionInterval column was called
# AnnualPeriodicity.
if ("AnnualPeriodicity" %in% names(Compadre)) {
  Compadre$ProjectionInterval <- Compadre$AnnualPeriodicity
}

comp_flag <- cdb_flag(Compadre, "check_NA_U")

comp_use <- subset(comp_flag, OrganismType == "Tree" &
  check_NA_U == FALSE &
  ProjectionInterval == 1)

## -----------------------------------------------------------------------------
CompadreData(comp_use)[, c(
  "SpeciesAccepted", "MatrixPopulation",
  "MatrixTreatment"
)]

## -----------------------------------------------------------------------------
# add column ID-ing matrices with same MatrixClassAuthor vector
comp_use$stage_id <- cdb_id_stages(comp_use)

# collapse database to single matrix per species * MatrixClassAuthor
comp_collapse <- cdb_collapse(comp_use, "stage_id")

# check species/populations again
CompadreData(comp_collapse)[, c(
  "SpeciesAccepted", "MatrixPopulation",
  "MatrixTreatment"
)]

## -----------------------------------------------------------------------------
MatrixClassOrganized(comp_collapse)

## -----------------------------------------------------------------------------
comp_collapse$start_life <- mpm_first_active(comp_collapse)

## ---- fig.width = 6, fig.height = 4-------------------------------------------
lx_list <- lapply(seq_len(nrow(comp_collapse)),
  function(x, comp_collapse) {
    U <- matU(comp_collapse$mat[[x]])

    rownames(U) <- colnames(U) # ensure row and col names are present

    mpm_to_lx(
      matU = U,
      start = comp_collapse$start_life[x],
      xmax = 40
    )
  },
  comp_collapse = comp_collapse
)

lx_array <- do.call(cbind, lx_list)

matplot(lx_array,
  type = "l", lty = 1, log = "y", ylim = c(0.0001, 1),
  lwd = 1.5, xlab = "Age (years)", ylab = "lx"
)

