# install.packages("BiocManager")
# BiocManager::install(c("graph", "RBGL", "Rgraphviz"))
# install.packages("pcalg")
library(pcalg)

set.seed(123)
p <- 10; w = 0.2
DAG <- randomDAG(p, w) ## true DAG
CPDAG <- dag2cpdag(DAG) ## true CPDAG
covTrue <- trueCov(DAG) ## true covariance matrix

par(mfrow = c(1,2))
Rgraphviz::plot(DAG)
Rgraphviz::plot(CPDAG)

  ## simulate Gaussian data from the true DAG
n <- 10000
X <- rmvDAG(n, DAG)

  ## estimate CPDAG and PDAG
suffStat <- list(C = cor(X), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, p=p, alpha = 0.01, conservative = TRUE)

  ## plot the true and estimated graphs
par(mfrow = c(1,3))
Rgraphviz::plot(DAG, main = "True DAG")
Rgraphviz::plot(CPDAG, main = "True CPDAG")
Rgraphviz::plot(pc.fit@graph, main = "Estimated CPDAG")

  ## Causal effect estimation using ida() function in pcalg
ida.out <- ida(3,10, crossprod(X), CPDAG, method = "local", type = "cpdag")
ida.out


# A different implementation of the IDA algorithm -------------------------

  ## Function that, for a given CPDAG, identifies all the possible distinct parent sets of node j
distinctPAfromCPDAG <- function(CPDAG, j) {
  CPDAGskel <- (CPDAG + t(CPDAG) != 0)*1
  pa <- setdiff(which(CPDAG[,j] != 0), which(CPDAG[j,] != 0))
  sib <- setdiff(which(CPDAG[,j] != 0), pa)
  nsib <- length(sib)
  if (nsib != 0) {
    if (nsib == 1) {
      pa.t <- list(integer(0), sib)
    } else {
      pa.t <- unlist(lapply(0:nsib, function(j) combn(sib, j, simplify = FALSE)), recursive = F)
    }
    ch.t <- lapply(1:length(pa.t), function(j) setdiff(sib, pa.t[[j]]))
    locallyvalidyn <- sapply(1:length(pa.t), function(i) {
      !pcalg:::has.new.coll(t(CPDAG), CPDAGskel, j, pa, pa.t[[i]], ch.t[[i]])
    })
    Z <- sapply(pa.t, function(i) c(pa, i))[locallyvalidyn]
  } else {
    Z <- pa
  }
  return(Z)
}

distinctPAfromCPDAG(as(CPDAG, "matrix"), 3)

  ## An implementation of the IDA algorithm that also returns the confidence intervals associated with each causal effect
ourIDA <- function(X, CPDAG, t, y, alpha) {
  validyn <- pcalg::isValidGraph(t(CPDAG), "cpdag")
  if (validyn == TRUE) {
    pas <- distinctPAfromCPDAG(CPDAG, t)
    ces <- vector("double", length(pas))
    confint <- matrix(0, nrow = length(pas), ncol = 2)
    if (length(pas) != 0) {
      for (j in 1:length(pas)) {
        if (!(y %in% pas[[j]])) {
          lmout <- lm(X[,y] ~ -1 + X[,c(t,pas[[j]])])
          ces[j] <- coef(lmout)[1]
          confint[j,] <- confint(lmout, level = alpha)[1,]
        }
      }
    } else {
      print("Input CPDAG is not valid")
      lmout <- lm(X[,y] ~ -1 + X[,t])
      ces <- coef(lmout)[1]
      confint <- confint(lmout, level = alpha)[1,]
    }
  } else {
    print("CPDAG not valid")
    pas <- integer(0)
    ces <- 0
    confint <- c(0,0)
  }
  return(list(ValidCPDAG = validyn, adj = pas, ce = ces, confint = confint))
}

ourIDA.out <- ourIDA(X, as(CPDAG, "matrix"), 3, 10, 0.95)
ourIDA.out$adj
ourIDA.out$ce
ourIDA.out$confint
