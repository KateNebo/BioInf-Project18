library("Rcpp")
library("gtools")
sourceCpp("scoreR.cpp")
source("se.R")
set.seed(42)

id <- 6 # id of peptide
df <- read.csv("data/Selected.csv", head = TRUE)[1:21,]
idx <- which(df$ID == id)
print(idx)

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
MetropolisUpdate_samc <- function(mass, score,tetta, beta)
{
      new.mass <- modify.mass(mass)
      new.score <- min(get.score(new.mass), MAX_SCORE)
      
      r <- exp((tetta[score+1] - tetta[new.score+1]) + beta*(score - new.score))
      alpha <- min(1,r)
      if (runif(1) < alpha) {
            score <- new.score
            mass <- new.mass
      }
      list(mass = mass, score = score)
}

update_tetta <- function(tetta, pi, gamma, l, s)
{
      for (i in 0:l)
      {
            tetta[i+1] = tetta[i+1] + gamma*((s == i) - pi[i+1])
      }
      
      return(tetta)
}

tetta <- rep(0, MAX_SCORE+1)
pi <- rep(1/(MAX_SCORE+1), MAX_SCORE+1)
ksi <- 0.8
t_0 <- 50
gamma <- unlist(lapply(seq(0,N), function(x) t_0/max(t_0, x^ksi)))
beta <- 0

# Sampling
s <- numeric(N)
mass <- as.numeric(rdirichlet(1, rep(1, N_MASS))*TOTAL_MASS)
score <- get.score(mass)
score <- min(score, MAX_SCORE)
s[1] <- score
for(j in 2:N)
{
      res <- MetropolisUpdate_samc(mass, s[j-1], tetta,beta)
      mass <- res$mass
      s[j] <- res$score
      
      tetta <- update_tetta(tetta, pi, gamma[j], MAX_SCORE, s[j])
}


# Conf intervals
const2 <- 1/N * sum(exp(tetta[s+1])*(pi[s+1]))
conf.samc2 <- function(x, beta, const, tetta)
{
      tmp <- se.bm(t(x), g <- function(x) (exp(tetta[x+1])*1/const*(pi[x+1])*(x >= MAX_SCORE)))
      se <- tmp$se.mean
      m <-tmp$mu
      z <- qnorm((1.95) / 2)
      return(list(m, m - z*se, m + z*se, se))
      
}



conf.int.samc2 <- conf.samc2(s, beta, const = const2, tetta = tetta)
res.samc <- unlist(conf.int.samc2)
samc.est = res.samc[1]
samc.se = res.samc[4]
