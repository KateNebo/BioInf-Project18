library("Rcpp")
library("gtools")
sourceCpp("scoreR.cpp")
source("se.R")
set.seed(42)


df <- read.csv("data/Selected.csv", head = TRUE)[1:21,]
id <- 6 # id of peptide
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
MetropolisUpdate <- function(mass, score, weight, Smin = 0)
{
      new.mass <- modify.mass(mass)
      new.score <- min(get.score(new.mass), MAX_SCORE)
      
      old.weight <- weight[score - Smin + 1]
      new.weight <- weight[new.score - Smin + 1]
      
      alpha <- min(1, exp(new.weight-old.weight))
      if (runif(1) < alpha) {
            score <- new.score
            mass <- new.mass
      }
      list(mass = mass, score = score)
}

isflat <- function(H, Smax, Smin)
{
      H1 <- H[Smin:Smax]
      m <- mean(H1)
      all(H1>20) & all(H1 > 0.75*m) & all(H1 < 1.25*m)
}

# Wanng-Landau algorithm
WangLandau <- function(Smin, Smax, w, phi =  0.1, phiFinal = 0.0002, thr = 10^5)
{
      mass <- as.numeric(rdirichlet(1, rep(1, N_MASS))*TOTAL_MASS)
      score <- get.score(mass)
      score <- min(score, MAX_SCORE)
      
      while(phi > phiFinal)
      {
            H <- rep(0, Smax)
            step = 0
            while(!isflat(H, Smax, Smin) & step < thr)
            {
                  res <- MetropolisUpdate(mass, score, w)
                  mass <- res$mass
                  score <- min(res$score, MAX_SCORE)
                  idx <- res$score - Smin + 1
                  H[idx] <- H[idx] + 1
                  w[idx] <- w[idx]-phi
                  step = step + 1
                  if (step %% 1000 == 0)
                  {
                        print(w)
                        print(H)
                  }
            }
            phi <- phi/2
            w <- w - max(w)
      }
      return(w)
}   


#Simple Monte Carlo
pval.est <- function(N) {
      res <- rdirichlet(N, rep(1, N_MASS))
      v <- apply(res, 1, function(x) get.score(x*TOTAL_MASS))
      v
}

v <- pval.est(N)
est <- length(v[v >= MAX_SCORE])/N


# Wang-Landau
s.min = 0
s.max = MAX_SCORE
w <-  rep(0, length.out = s.max - s.min + 1)

weights <- WangLandau(Smin = 0, Smax = MAX_SCORE, w)

s <- numeric(N)
mass <- as.numeric(rdirichlet(1, rep(1, N_MASS))*TOTAL_MASS)
score <- get.score(mass)
score <- min(score, MAX_SCORE)
s[1] <- score
for(j in 2:N)
{
      
      res <- MetropolisUpdate(mass, s[j-1], weights)
      mass <- res$mass
      s[j] <- res$score
}

resWL <- s
hist <- table(resWL)
const.wl <- sum(exp(-weights[resWL+1]))/N

#Conf intervals
conf.w <- function(x, w, const)
{
      tmp <- se.bm(t(x), g <- function(x) (x >= MAX_SCORE) * exp(-w[x +1]) / const)
      se <- tmp$se.mean
      m <-tmp$mu
      z <- qnorm((1.95) / 2)
      res <- c(m, m - z*se, m + z*se, se)
      res
}



res.wl <- conf.w(resWL, weights, const.wl)

wl.est = res.wl[1]
wl.se = res.wl[4]
mc.est = est
mc.se = sqrt(est*(1 - est)/N)

