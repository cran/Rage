## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setupDarwin, include=FALSE, eval = Sys.info()[["sysname"]] == "Darwin"----
#The following line seems to be required by pkgdown::build_site() on my machine, but causes build to break with R-CMD-CHECK on GH
knitr::opts_chunk$set(dev = "png", dev.args = list(type = "cairo-png"))

## ---- message=FALSE-----------------------------------------------------------
library(Rcompadre)
library(popdemo)
library(Rage)
library(ggtern)

## ----load example COMPADRE data, eval=TRUE------------------------------------
data(Compadre)
Compadre <- cdb_flag(Compadre)
Compadre <- subset(Compadre, MatrixSplit == "Divided" & check_ergodic == TRUE)

## ----elas, eval=TRUE----------------------------------------------------------
popdemo::elas(matA(Compadre)[[1]])

## -----------------------------------------------------------------------------
mat_a <- matA(Compadre)[[1]]

perturb_matrix(mat_a, type = "elasticity")

#Rage::matrixElementPerturbation(matU = matU(Compadre)[[1]],    This function used to be able to accommodate partioning
#                                matF = matF(Compadre)[[1]],    into U F and C so we coudl partition this... but it's now gone
#                                matC = matC(Compadre)[[1]])    in perturb_matrix.... why?

## -----------------------------------------------------------------------------
Amats <- matA(Compadre)
Umats <- matU(Compadre)
Fmats <- matF(Compadre)
Cmats <- matC(Compadre)

output <- data.frame(S=rep(NA,length(Umats)),G=NA,R=NA,lam=NA)

for(i in 1:length(Umats)){
  temp <- perturb_vr(Umats[[i]],   #The same issue as above
                     Fmats[[i]],
                     Cmats[[i]],
                     type = "elasticity")
  
  output$S[i] <- temp$survival
  output$G[i] <- temp$growth + temp$shrinkage
  output$R[i] <- temp$fecundity + temp$clonality
  
  #Calculate growth rate
  output$lam[i] <- popdemo::eigs(Amats[[i]], "lambda")
}

## -----------------------------------------------------------------------------
head(output)

## -----------------------------------------------------------------------------
output[,1:3] <- t(apply(output[,1:3], 1, function(x) x/sum(x)))

## -----------------------------------------------------------------------------
B<-ggtern::ggtern(data = output,
                  aes(x = R,
                      y = G,
                      z = S,
                      colour = lam))  +  
  geom_point() +    
  scale_color_viridis_c()+
  theme_showarrows() +
  theme_clockwise() 

B

## -----------------------------------------------------------------------------

data(mpm1)

perturb_vr(mpm1$matU,
           mpm1$matF,
           mpm1$matF,
           type = "elasticity")

## -----------------------------------------------------------------------------
Amats <- matA(Compadre)
Umats <- matU(Compadre)
Fmats <- matF(Compadre)
Cmats <- matC(Compadre)

output <- data.frame(survival=rep(NA,length(Umats)),growthShrinkage=NA,reproduction=NA,lam=NA)
for(i in 1:length(Umats)){
  temp <- perturb_vr(Umats[[i]],
                     Fmats[[i]],
                     Cmats[[i]],
                     type = "elasticity")
  
output$survival[i] <- temp$survival
output$growthShrinkage[i] <- temp$growth + temp$shrinkage
output$reproduction[i] <- temp$fecundity + temp$clonality

#Calculate growth rate
output$lam[i] <- popdemo::eigs(Amats[[i]], "lambda")
}

## -----------------------------------------------------------------------------
head(output)

## -----------------------------------------------------------------------------
output[,1:3] <- t(apply(output[,1:3], 1, function(x) x/sum(x)))

## -----------------------------------------------------------------------------

B<-ggtern::ggtern(data = output,
                  aes(x = reproduction,
                      y = growthShrinkage,
                      z = survival,
                      colour = lam))  +  
  geom_point() +  
  scale_color_viridis_c()+
  theme_showarrows() +
  theme_clockwise() 

B

