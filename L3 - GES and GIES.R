# install.packages("BiocManager")
# BiocManager::install(c("graph", "RBGL", "Rgraphviz"))
# install.packages("pcalg")
# install.packages("BCDAG")
library(pcalg)


# Generate DAG, parameters and data ---------------------------------------

p <- 10; w <- 0.15
set.seed(1)
DAG1 <- BCDAG::rDAG(p,w)
L1 <- diag(1,p) - runif(p^2, 1, 2)*DAG1*sample(c(-1,1), p^2, replace = TRUE)
Sigma1 <- crossprod(solve(L1))

n <- 10^3
set.seed(1)
X1 <- mvtnorm::rmvnorm(n, sigma = Sigma1)


# Use GES -----------------------------------------------------------------

# ?pcalg::ges
# ?pcalg::`Score-class`

score <- new("GaussL0penObsScore", data = X1) #Using BIC
ges_out <- ges(score)

par(mfrow = c(1,3))
Rgraphviz::plot(BCDAG::as_graphNEL(DAG1), main = "True DAG")
Rgraphviz::plot(BCDAG::as_graphNEL(dag2essgraph(DAG1)), main = "True CPDAG")
Rgraphviz::plot(ges_out$essgraph, main = "Estimated via GES")


# Generate interventional data --------------------------------------------

t <- 4
L1_int <- L1;
L1_int[,t] <- 0; L1_int[t,t] <- 1
Sigma1_int <- crossprod(solve(L1_int))

set.seed(1)
X1_int <- mvtnorm::rmvnorm(n/2, sigma = Sigma1_int)

# Use GIES ----------------------------------------------------------------

?pcalg::gies
# ?pcalg::`Score-class`

targets <- list(integer(0), t)
target.index <- c(rep(1,n), rep(2,n/2))
X <- rbind(X1, X1_int)
score_int <- new("GaussL0penIntScore", data = X, targets = targets, target.index = target.index)
gies_out <- gies(score_int)

par(mfrow = c(1,3))
Rgraphviz::plot(BCDAG::as_graphNEL(DAG1), main = "True DAG")
Rgraphviz::plot(BCDAG::as_graphNEL(dag2essgraph(DAG1, targets = list(integer(0), 4))), main = "True I-Ess")
Rgraphviz::plot(gies_out$essgraph, main = "Estimated via GIES")

