## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(Rage)

## -----------------------------------------------------------------------------
## survival probabilities
s1 <- 0.20 # seed
s2 <- 0.40 # small
s3 <- 0.50 # large
s4 <- 0.30 # dormant

## matU vital rates conditional on survival (growth, shrinkage, stasis, etc.)
g12 <- 0.15 # seed to small      (growth)
g13 <- 0.05 # seed to large      (growth)
g11 <- 1 - g12 - g13 # seed to seed       (stasis)
g23 <- 0.45 # small to large     (growth)
g24 <- 0.10 # small to dormant   (enter dormancy)
g22 <- 1 - g23 - g24 # small to small     (stasis)
g34 <- 0.22 # large to dormant   (enter dormancy)
g33 <- 1 - g34 # large to large     (stasis)
g43 <- 0.50 # dormant to large   (exit dormancy)
g44 <- 1 - g43 # dormant to dormant (stasis)

## matF vital rates conditional on survival (i.e. fecundity)
f2 <- 0.4 # small
f3 <- 1.1 # large

## -----------------------------------------------------------------------------
# growth/survival component
matU <- rbind(
  c(s1 * g11, 0, 0, 0),
  c(s1 * g12, s2 * g22, 0, 0),
  c(s1 * g13, s2 * g23, s3 * g33, s4 * g43),
  c(0, s2 * g24, s3 * g34, s4 * g44)
)

# sexual reproduction component
matF <- rbind(
  c(0, s2 * f2, s3 * f3, 0),
  c(0, 0, 0, 0),
  c(0, 0, 0, 0),
  c(0, 0, 0, 0)
)

## -----------------------------------------------------------------------------
(surv <- vr_vec_survival(matU))
# equivalent to...
(surv <- colSums(matU))

## -----------------------------------------------------------------------------
# matU vital rates conditional on survival
vr_mat_U(matU)

# matF vital rates conditional on survival
vr_mat_R(matU, matF)

## -----------------------------------------------------------------------------
vr_vec_growth(matU)

## -----------------------------------------------------------------------------
vr_vec_growth(matU, exclude_row = 4)

## -----------------------------------------------------------------------------
vr_vec_shrinkage(matU)

## -----------------------------------------------------------------------------
vr_vec_shrinkage(matU, exclude_col = 4)

## -----------------------------------------------------------------------------
vr_vec_stasis(matU)

## -----------------------------------------------------------------------------
vr_vec_dorm_enter(matU, dorm_stages = 4)
vr_vec_dorm_exit(matU, dorm_stages = 4)

## -----------------------------------------------------------------------------
vr_vec_reproduction(matU, matF)

## -----------------------------------------------------------------------------
vr_growth(matU, exclude_row = 4)

# equivalent to
(vec_growth <- vr_vec_growth(matU, exclude_row = 4))
mean(vec_growth, na.rm = TRUE)

## -----------------------------------------------------------------------------
# calculate the stable distribution using popdemo::eigs
library(popdemo)
matA <- matU + matF
w <- popdemo::eigs(matA, what = "ss")

# calculate mean vital rate of growth weighted by the stable distribution
vr_growth(matU, exclude_row = 4, weights_col = w)

# equivalent to
pos <- !is.na(vec_growth) # possible transitions
sum(vec_growth[pos] * w[pos]) / sum(w[pos]) # weighted average

## -----------------------------------------------------------------------------
vr_vec_shrinkage(matU)
vr_shrinkage(matU)

## -----------------------------------------------------------------------------
(pos <- matU > 0) # possible transitions based on matU
pos[2, 3] <- TRUE # set medium-to-small shrinkage tr to possible

vr_vec_shrinkage(matU, posU = pos)
vr_shrinkage(matU, posU = pos)

