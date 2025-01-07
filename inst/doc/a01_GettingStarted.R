## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE,fig.align='center',fig.height=2,fig.width=6-------------------
# hidden code to produce figures
library(DiagrammeR)
matA <- rbind(
  c(0.0, 0.0, 3.2),
  c(0.5, 0.3, 0.8),
  c(0.0, 0.4, 0.9)
)
stages <- c("seedling", "rosette", "flowering")
title <- NULL
graph <- expand.grid(to = stages, from = stages)
graph$trans <- round(c(matA), 3)
graph <- graph[graph$trans > 0, ]
nodes <- paste(paste0("'", stages, "'"), collapse = "; ")
graph$min_len <- (as.numeric(graph$to) - as.numeric(graph$from)) * 3
graph$col <- c(
  "PaleGreen4", "PaleGreen4", "PaleGreen4", "Goldenrod1",
  "MediumOrchid4", "PaleGreen4"
)
edges <- paste0("'", graph$from, "'", " -> ", "'", graph$to, "'",
  "[minlen=", graph$min_len,
  ",fontsize=", 10,
  ",color=", graph$col,
  ",xlabel=", paste("\"", graph$trans),
  "\"]\n",
  collapse = ""
)
grViz(
  paste(
    "
digraph {
  {
    graph[overlap=false];
    rank=same;
    node [shape=", "egg", ", fontsize=", 12, "];",
    nodes, "
  }",
    "ordering=out
  x [style=invis]
  x -> {", nodes, "} [style=invis]", edges,
    "labelloc=\"t\";
  label=\"", title, "\"
}"
  )
)

## ----echo=FALSE,fig.align='center',fig.height=4,fig.width=4-------------------
library(ggplot2)
ggdat <- merge(graph,
  expand.grid(to = stages, from = stages),
  by = c("to", "from"),
  all.x = TRUE, all.y = TRUE
)
ggdat$trans[is.na(ggdat$trans)] <- 0
ggdat$col[is.na(ggdat$col)] <- "transparent"
ggdat$to <- factor(ggdat$to, levels = c("flowering", "rosette", "seedling"))
ggdat$from <- factor(ggdat$from, levels = c("seedling", "rosette", "flowering"))
ggplot(ggdat, aes(x = from, y = to, label = trans)) +
  geom_tile(color = "black", fill = "white", linewidth = 0.25, show.legend = FALSE) +
  geom_text(size = 6) +
  scale_x_discrete(position = "top") +
  labs(x = "current life stage", y = "life stage at time t+1") +
  coord_equal(expand = FALSE) +
  theme_bw(base_size = 18) +
  theme(panel.border = element_blank())

## ----echo=FALSE,fig.align='center',fig.height=4,fig.width=9,out.width='100%'----
blankdat <- expand.grid(to = stages, from = stages, trans = 0)
blankdat$to <- factor(blankdat$to, levels = c(
  "flowering", "rosette",
  "seedling"
))
blankdat$from <- factor(blankdat$from,
  levels = c("seedling", "rosette", "flowering")
)
ggdat$col <- factor(ggdat$col,
  levels = c(
    "PaleGreen4", "Goldenrod1",
    "MediumOrchid4", "transparent"
  ),
  labels = c("U", "F", "C", "t")
)
ggplot(
  ggdat[ggdat$col != "t", ],
  aes(x = from, y = to, fill = col, label = trans)
) +
  geom_tile(
    data = blankdat,
    aes(fill = NULL), fill = "white", color = "black", linewidth = 0.25
  ) +
  geom_text(data = blankdat, aes(fill = NULL), size = 6) +
  geom_tile(color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(size = 6) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c(
    "F" = "goldenrod1",
    "C" = "mediumorchid4",
    "U" = "palegreen4",
    "t" = "white"
  )) +
  labs(x = "current life stage", y = "life stage at time t+1") +
  coord_equal(expand = FALSE) +
  facet_wrap(~col, nrow = 1) +
  theme_bw(base_size = 18) +
  theme(
    panel.border = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside"
  )

## ----message=FALSE------------------------------------------------------------
library(popdemo)

# define the transition matrix, A
A <- rbind(
  c(0.0, 0.0, 3.2),
  c(0.5, 0.3, 0.8),
  c(0.0, 0.4, 0.9)
)

# lambda: equilibrium per-capita population growth rate
popdemo::eigs(A = A, what = "lambda")

# w: stable stage distribution (relative frequencies)
popdemo::eigs(A = A, what = "ss")

## -----------------------------------------------------------------------------
library(Rage) # load Rage
data(mpm1) # load data object 'mpm1'
mpm1 # display the contents

## ----echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----------------
tabl <- "
| Function category                       | Stand-alone vignette             |
|-----------------------------------------|----------------------------------|
| 1. [Vital rates](#vitalrates)  | [VitalRates](a02_VitalRates.html)     |
| 2. [Life tables](#lifetable) | [AgeFromStage](a04_AgeFromStage.html) |
| 3. [Perturbation analysis](#perturb)    | n/a |
| 4. [Deriving life history traits](#lifehist) | [LifeHistoryTraits](a03_LifeHistoryTraits.html) |
| 5. [Transformation of matrices](#maniptransform) | n/a                                   |
"
cat(tabl)

## -----------------------------------------------------------------------------
vr_vec_survival(matU = mpm1$matU)
vr_vec_stasis(matU = mpm1$matU)

## -----------------------------------------------------------------------------
# product of Pr(survival) and Pr(stasis) yields Pr(stasis|survived)
vr_vec_survival(matU = mpm1$matU) * vr_vec_stasis(matU = mpm1$matU)
diag(mpm1$matU) # equivalent to the diagonal of U matrix

## -----------------------------------------------------------------------------
vr_survival(matU = mpm1$matU, exclude_col = 1) # exclude 'seed' stage
mean(vr_vec_survival(mpm1$matU)[-1]) # equivalent to the mean without 'seed'

## -----------------------------------------------------------------------------
lt <- mpm_to_table(matU = mpm1$matU, matF = mpm1$matF) # full life table
lt
lx <- mpm_to_lx(matU = mpm1$matU) # survivorship to start of each age class
lx

## -----------------------------------------------------------------------------
lx_to_px(lx = lx) # survivorship to survival probability
lx_to_hx(lx = lx) # survivorship to mortality hazard

## ----warning=FALSE, message=FALSE, fig.align='center', fig.height=5, fig.width=6----
# project a germinated cohort through the U matrix
cohort <- popdemo::project(A = mpm1$matU, vector = c(0, 1, 0, 0, 0), time = 10)
popStructure <- vec(cohort) / rowSums(vec(cohort))

matplot(popStructure,
  type = "l", xlab = "time",
  ylab = "proportion in stage class"
)

## -----------------------------------------------------------------------------
# calculate time to QSD from the U matrix of an MPM
(q <- qsd_converge(mat = mpm1$matU, start = "small"))

# subset the life table rows to ages prior to the QSD
lt_preQSD <- lt[1:q, ]

# plot mortality trajectory from the life table subset (blue),
# showing plateau effect if the trajectory (grey) was allowed to continue to the
# QSD (dashed vertical line) and beyond
plot(qx ~ x,
  data = lt, type = "l", col = "darkgrey", ylim = c(0, 1),
  xlab = "age"
)
lines(qx ~ x, data = lt_preQSD, type = "l", col = "blue", lwd = 4)
abline(v = q, lty = "dashed")

## -----------------------------------------------------------------------------
# construct the transition matrix A = U + F (+ C when present)
mpm1$matA <- with(mpm1, matU + matF)

# sensitivity of lambda to...
# ...matrix element perturbations
perturb_matrix(
  matA = mpm1$matA,
  type = "sensitivity", demog_stat = "lambda"
)
# ...vital rate perturbations
perturb_vr(
  matU = mpm1$matU, matF = mpm1$matF,
  type = "sensitivity", demog_stat = "lambda"
)
# ...transition type perturbations
perturb_trans(
  matU = mpm1$matU, matF = mpm1$matF,
  type = "sensitivity", demog_stat = "lambda"
)

## -----------------------------------------------------------------------------
# post-germination time steps until post-germination survivorship falls below 5%
longevity(matU = mpm1$matU, start = "small", lx_crit = 0.05)
# expected lifetime production of 'small' offspring by a 'small' individual
net_repro_rate(
  matU = mpm1$matU, matR = mpm1$matF, start = "seed",
  method = "start"
)

## -----------------------------------------------------------------------------
# derive post-germination survivorship trajectory from U matrix
lx <- mpm_to_lx(matU = mpm1$matU, start = "small")

## -----------------------------------------------------------------------------
# collapse 'small', 'medium', and 'large' stages into single stage class
col1 <- mpm_collapse(
  matU = mpm1$matU, matF = mpm1$matF,
  collapse = list(1, 2:4, 5)
)
col1$matA

## -----------------------------------------------------------------------------
# automated stage naming
(col1_auto <- name_stages(mat = col1, prefix = "class_"))

# overwrite with custom stages
(col1_cust <- name_stages(
  mat = col1, names = c("seed", "active", "dormant"),
  prefix = NULL
))

## -----------------------------------------------------------------------------
# compare population growth rate of original and collapsed MPM (preserved)
popdemo::eigs(A = mpm1$matA, what = "lambda")
popdemo::eigs(A = col1_cust$matA, what = "lambda")

# compare net reproductive rate of original and collapsed MPM (not preserved)
net_repro_rate(matU = mpm1$matU, matR = mpm1$matF)
net_repro_rate(matU = col1_cust$matU, matR = col1_cust$matF)

