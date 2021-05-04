## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setupDarwin, include=FALSE, eval = Sys.info()[["sysname"]] == "Darwin"----
#The following line seems to be required by pkgdown::build_site() on my machine, but causes build to break with R-CMD-CHECK on GH
knitr::opts_chunk$set(dev = "png", dev.args = list(type = "cairo-png"))

## -----------------------------------------------------------------------------
library(Rage)  # load Rage
data(mpm1)     # load data object 'mpm1'
mpm1           # display the contents

## -----------------------------------------------------------------------------
life_expect_mean(matU = mpm1$matU, start = 1)  # mean life expectancy from "seed" stage
life_expect_var(matU = mpm1$matU, start = 1)   # variance in life expectancy from "seed" stage

## -----------------------------------------------------------------------------
life_expect_mean(matU = mpm1$matU, start = 2)  # mean life expectancy from "small" stage

## -----------------------------------------------------------------------------
life_expect_mean(matU = mpm1$matU, start = c(0,0.4,0.6,0,0))    

## -----------------------------------------------------------------------------
longevity(matU = mpm1$matU, start = 2, lx_crit = 0.05)  

## -----------------------------------------------------------------------------
longval <- c()
startvec <- c(1:dim(mpm1$matU)[1])   # vector of starting stages

for(i in c(startvec)) {
  longval[i] <- longevity(matU = mpm1$matU, start = startvec[i], lx_crit = 0.05)  
}
  
plot(longval, type = "l", xlab = "Starting stage", ylab = "Longevity to 5% survivorship")
longval

## -----------------------------------------------------------------------------
net_repro_rate(matU = mpm1$matU, matR = mpm1$matF)            

## -----------------------------------------------------------------------------
gen_time(matU = mpm1$matU, matR = mpm1$matF)    

## -----------------------------------------------------------------------------
mature_age(matU = mpm1$matU, matR = mpm1$matF, start = 2)     

## -----------------------------------------------------------------------------
mature_prob(matU = mpm1$matU, matR = mpm1$matF, start = 2)   

## -----------------------------------------------------------------------------
mpm1$matF # We see that the "medium" and "large" stages are reproductive

maturedist <- mature_distrib(matU = mpm1$matU, start = 1L, repro_stages = c(FALSE,FALSE,TRUE,TRUE,FALSE)) 
maturedist

## -----------------------------------------------------------------------------
life_expect_mean(matU = mpm1$matU, start = c(maturedist))   # mean life expectancy from maturity
life_expect_mean(matU = mpm1$matU, start = 2)       # mean life expectancy from "small" stage

## -----------------------------------------------------------------------------
lx <- mpm_to_lx(matU = mpm1$matU, start = "small")
px <- mpm_to_px(matU = mpm1$matU, start = "small")
hx <- mpm_to_hx(matU = mpm1$matU, start = "small")
mx <- mpm_to_mx(matU = mpm1$matU, matR = mpm1$matF, start = "small")

## -----------------------------------------------------------------------------
lx_seed <- mpm_to_lx(matU = mpm1$matU, start = "seed")

plot(lx, xlab="Survival time (years)", ylab="Survivorship", type="s", col="black")
  lines(lx_seed, type="s", col="orange")
  legend("topright", inset=c(0.05,0.05), c("From seed stage", "From small stage"), 
       lty = c(1,1), col = c("orange","black"))

## -----------------------------------------------------------------------------
entropy_d(lx, mx)   # Demetrius' entropy
entropy_k(lx)       # Keyfitz' entropy

## -----------------------------------------------------------------------------
shape_surv(lx)      # shape of survival/mortality trajectory
shape_rep(mx)       # shape of fecundity trajectory
