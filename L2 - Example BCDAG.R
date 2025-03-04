# install.packages("BiocManager")
# BiocManager::install(c("graph", "Rgraphviz", "RBGL"))
# install.packages("BCDAG")
library(BCDAG)


###################
## Gaussian data ##
###################

##################################
## DAG and parameter simulation ##
##################################

p = 8

set.seed(1)
dag = rDAG(p, w = 0.2)

Rgraphviz::plot(as(dag, "graphNEL"))

set.seed(123)
DL = rDAGWishart(n = 1, DAG = dag, a = 10, U = diag(1, p))

D = DL$D
L = DL$L

#####################
## Data generation ##
#####################

Sigma = solve(t(L))%*%D%*%solve(L)

set.seed(123)
X = mvtnorm::rmvnorm(n = 200, mean = rep(0, 8), sigma = Sigma)

## Setting prior hyperparameters

p = ncol(X)
n = nrow(X)

a = p
U = diag(1,p)

#####################################
## Algorithm 1 - Collapsed sampler ##
#####################################

?learn_DAG

## Run learn_DAG with setting collapse = TRUE
S <- p*1000; burn = S/2
out_alg_1 = learn_DAG(S = S, burn = burn, data = X, a = a, U = U, w = 0.2,
                      fast = TRUE, collapse = TRUE)

## Available methods for posterior summaries

print(out_alg_1)
summary(out_alg_1)
plot(out_alg_1)

#################################
## Algorithm 2 - Joint sampler ##
#################################

## Run learn_DAG with setting collapse = FALSE

out_alg_2 = learn_DAG(S = S, burn = burn, data = X, a = a, U = U, w = 0.2,
                      fast = TRUE, collapse = FALSE)

print(out_alg_2)
summary(out_alg_2)
plot(out_alg_2)

get_diagnostics(out_alg_2)

## Recover (get) the posterior of causal effect parameter of do(X_4 = x) on Y = X_1

causal_4 = get_causaleffect(out_alg_2, targets = 4, response = 1)

## Available methods for posterior summaries of causal effects

summary(causal_4)
plot(causal_4)




