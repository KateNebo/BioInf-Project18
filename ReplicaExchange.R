library("Rcpp")
library("gtools")
sourceCpp("scoreR.cpp")
source("se.R")
set.seed(42)

id <- 6 # id of peptide
df <- read.csv("data/Selected.csv", head = TRUE)[1:21,]
idx <- which(df$ID == id)

mat <- as.matrix(read.table(paste0("data/mgf_input/matrix/matrix", df$LENGTH[idx], ".txt")))
rule <- as.matrix(read.table(paste0("data/mgf_input/rule/rule", df$LENGTH[idx], ".txt")))
exp_spectrum <- as.numeric(read.table(paste0("data/mgf_input/spectrum/spectrum", id, ".txt")))

MASS_PROTON = 1.00728
N_MASS = ncol(mat)
TOTAL_MASS = df$MASS[idx]
MAX_SCORE = df$SCORE.MDDPR.[idx]
N = 10^7

exp_spectrum <- sort(exp_spectrum)

modify.mass <- function(mass)
      update_mass(mass, rule)

get.score <- function(mass)
      score_peak(exp_spectrum, mat %*% mass, TOTAL_MASS, MASS_PROTON, FALSE)

# Metropolis-Hastings algorithm
MetropolisUpdate_re <- function(mass, score, beta)
{
      new.mass <- modify.mass(mass)
      new.score <- min(get.score(new.mass), MAX_SCORE)
      alpha <- min(1,exp(beta*(new.score - score)))
      
      if (runif(1) < alpha) {
            score <- new.score
            mass <- new.mass
      }
      list(mass = mass, score = score)
}

beta <- seq(0,1, 0.1)
step_swap <- 500
sampling <- T
num_swap <- 0
k <-length(beta)

# Sampling
if(sampling)
{
  s <- matrix(0, k, N)
  mass <- replicate(k, as.numeric(rdirichlet(1, rep(1, N_MASS))*TOTAL_MASS))
  s[,1] <- unlist(apply(mass,2, get.score))
  s[,1] <- unlist(lapply(s[,1], min, MAX_SCORE))
  for(j in 2:N)
  {
    for(i in 1:k)
    {
      res <- MetropolisUpdate_re(mass[,i], s[i,j-1], beta[i])
      mass[,i] <- res$mass
      s[i,j] <- res$score
    }
    
    if (j%%step_swap ==0)
    {
      k_swap = sample(1:k, 2)
      alpha_swap = min(1, exp((beta[k_swap[1]] - beta[k_swap[2]])* (s[k_swap[1], j] - s[k_swap[2], j])))
      if(runif(1)< alpha_swap)
      {
        tmp = mass[,k_swap[1]]
        mass[,k_swap[1]] = mass[,k_swap[2]]
        mass[,k_swap[2]] = tmp
        num_swap <- num_swap + 1
      }
    }
  }
  
}


# Conf intervals
conf.re <- function(x, beta)
{
  const <- 1/N * sum(exp(-beta*x))
  tmp <- se.bm(t(x), g <- function(x) (x >= MAX_SCORE)*((exp(-beta*x)*1/const)))
  se <- tmp$se.mean
  m <-tmp$mu
  z <- qnorm((1.95) / 2)
  return(list(m, m - z*se, m + z*se, se))
  
}

ci.re <- sapply(c(1:k), function(j) conf.re(s[j,],beta[j]))

res.df = cbind(beta, as.data.frame(t(ci.re)))

m = mean(unlist(res.df$V1))
se = sqrt(sum(unlist(res.df$V4)^2)/k^2)
