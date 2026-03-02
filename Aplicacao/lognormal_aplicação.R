options(scipen = 999)  # evita notação científica
set.seed(1)

library(survival)
library(nlme)   # para fdHess
library(MASS)   # para ginv

##-----------------------------------------------------------------
# Dados Infarto Agudo do Miocárdio
##-----------------------------------------------------------------
dados <- read.table("dados-infarto.txt")

# Tempo e status
y      <- dados$fim - dados$ini
status <- dados$status

# Covariáveis: idade (<65 / ≥65) e sexo
x1 	   <- factor(ifelse(dados$idade < 65, "(<65)", "(>=65)"))
x2 	   <- factor(dados$sexo)
dados      <- data.frame(dados, x1, x2)

#-----------------------------------------------------------------
# Matrizes de covariáveis
#-----------------------------------------------------------------
Xe_sigma <- model.matrix(~ 1, data = dados)         # para σ
Xe_mu    <- model.matrix(~ 1 + x1 + x2, data = dados)  # para μ
Xe_p0    <- model.matrix(~ 1 + x1 + x2, data = dados)  # para p₀
Xe_phi   <- model.matrix(~ 1, data = dados)         # para φ

#-----------------------------------------------------------------
# Função de verossimilhança: Lognormal Discreta + fração de cura
#-----------------------------------------------------------------
fvero <- function(theta) {
  n_sigma <- ncol(Xe_sigma)
  n_mu    <- ncol(Xe_mu)
  n_p0    <- ncol(Xe_p0)
  n_phi   <- ncol(Xe_phi)

  # índices dos parâmetros
  beta_sigma <- theta[1:n_sigma]
  beta_mu    <- theta[(n_sigma+1):(n_sigma+n_mu)]
  beta_p0    <- theta[(n_sigma+n_mu+1):(n_sigma+n_mu+n_p0)]
  beta_phi   <- theta[(n_sigma+n_mu+n_p0+1):(n_sigma+n_mu+n_p0+n_phi)]

  # transformações
  sigma <- exp(Xe_sigma %*% beta_sigma)
  mu    <- Xe_mu %*% beta_mu
  p0    <- 1/(1 + exp(-Xe_p0 %*% beta_p0))
  phi   <- exp(Xe_phi %*% beta_phi)
  eta   <- phi^(-1) * (p0^(-phi) - 1)

  # Função de sobrevivência Lognormal Discreta
  S_LND <- function(y, sigma, mu) {
    1 - pnorm((log(y) - mu) / sigma)
  }

  vF1   <- 1 - S_LND(y, sigma, mu)
  Spop1 <- (1 + phi * eta * vF1)^(-1/phi)

  vF2   <- 1 - S_LND(y + 1, sigma, mu)
  Spop2 <- (1 + phi * eta * vF2)^(-1/phi)

  fpop  <- Spop1 - Spop2

  # log-verossimilhança
  loglik <- sum(status * log(fpop) + (1 - status) * log(Spop1))
  return(loglik)
}

#-----------------------------------------------------------------
# Estimação
#-----------------------------------------------------------------
n_sigma <- ncol(Xe_sigma)
n_mu    <- ncol(Xe_mu)
n_p0    <- ncol(Xe_p0)
n_phi   <- ncol(Xe_phi)

theta0 <- rep(0, n_sigma + n_mu + n_p0 + n_phi)

mgg1 <- optim(theta0, fvero, method = "SANN",
              control = list(fnscale = -1, maxit = 10000))
theta0 <- mgg1$par

mgg <- optim(theta0, fvero, method = "BFGS",
             control = list(fnscale = -1, maxit = 1000))

#-----------------------------------------------------------------
# Parâmetros estimados
#-----------------------------------------------------------------
beta_sigma <- mgg$par[1:n_sigma]
beta_mu    <- mgg$par[(n_sigma+1):(n_sigma+n_mu)]
beta_p0    <- mgg$par[(n_sigma+n_mu+1):(n_sigma+n_mu+n_p0)]
beta_phi   <- mgg$par[(n_sigma+n_mu+n_p0+1):(n_sigma+n_mu+n_p0+n_phi)]

sigma <- exp(Xe_sigma %*% beta_sigma)
mu    <- Xe_mu %*% beta_mu
p0    <- 1/(1 + exp(-Xe_p0 %*% beta_p0))
phi   <- exp(Xe_phi %*% beta_phi)
eta   <- phi^(-1) * (p0^(-phi) - 1)

S_LND <- function(y, sigma, mu) {
  1 - pnorm((log(y) - mu) / sigma)
}

vF   <- 1 - S_LND(y + 1, sigma, mu)
Spop <- (1 + phi * eta * vF)^(-1/phi)
p0   <- (1 + phi * eta)^(-1/phi)

#-----------------------------------------------------------------
# Plot função de sobrevivência empírica vs modelo
#-----------------------------------------------------------------
mKM <- survfit(Surv(y, status) ~ x1 + x2, se.fit = FALSE)
plot(mKM, xlab = "Tempo (meses)", ylab = "Função de sobrevivência",
     cex.axis = 1.5, cex.lab = 1.5)
points(y, Spop, pch = 20, cex = 0.8)

#-----------------------------------------------------------------
# Erros padrão (pseudo-inversa da Hessiana)
#-----------------------------------------------------------------

Estimativa <- mgg$par
obsinf <- -fdHess(Estimativa, fvero)$Hessian
covmat <- ginv(obsinf)
setheta <- sqrt(pmax(0, diag(covmat)))
tvals <- Estimativa / setheta
pvals <- 2 * (1 - pnorm(abs(tvals)))

# nomes dos parâmetros
nomes <- c(paste0("b_sigma_", colnames(Xe_sigma)),
           paste0("b_mu_",    colnames(Xe_mu)),
           paste0("b_p0_",    colnames(Xe_p0)),
           paste0("b_phi_",   colnames(Xe_phi)))

mest <- cbind(
  Estimate = Estimativa,
  `s.e.` = setheta,
  `|t|` = abs(tvals),
  `p-value` = pvals
)
rownames(mest) <- nomes

print(mest, 6)

cormat <- cov2cor(covmat)
print(cormat, 3)

#-----------------------------------------------------------------
# AIC
#-----------------------------------------------------------------
2 * length(Estimativa) - 2 * fvero(Estimativa)

# Teste Log-Rank
teste_log_rank <- survdiff(Surv(y, status) ~ x1 + x2)
print(teste_log_rank)

