library(survival)
library(nlme)    # Para fdHess (Hessiana)
library(MASS)    # Para ginv (Pseudo-inversa)
library(stats)   # Para uniroot, pnorm, qnorm, optim, rnbinom
library(ggplot2) # Para gráficos
library(tidyr)   # Para pivot_longer (processamento de dados para gráficos)
library(extraDistr)
set.seed(1)

# =================================================================
# SECTION 0: DEFINIÇÕES E CONFIGURAÇÃO
# =================================================================

# Parâmetros verdadeiros (True Values)
PARAM_FIXOS_ADJ <- list(
    beta_sigma_int = 0.577998, 
    beta_mu_int = 2.753863,
    beta_mu_x1 = 1.937027,
    beta_mu_x2 = -0.887691,
    beta_p0_int = 1.688387,
    beta_p0_x1 = -1.206960,
    beta_p0_x2 = 0.420279,
    beta_phi_int = 2.297517 
)

THETA_TRUE <- unlist(PARAM_FIXOS_ADJ)
N_PARAM <- length(THETA_TRUE) 
names(THETA_TRUE) <- names(PARAM_FIXOS_ADJ)

# Variáveis do estudo de simulação
SAMPLE_SIZES <- c(200)
N_MONTE_CARLO <- 250
P_SUSC_TARGET <- 0.05 


# =================================================================
# SECTION 1: FUNÇÕES LND DISCRETA
# =================================================================


# Geração de Tempo Aleatório LND Discreto
rlnd_disc <- function(n, mu, sigma) {
    U <- runif(n)
    Z <- qnorm(U)
    W <- ceiling(1+exp(mu + sigma * Z))
    return(pmax(0, W))
}

# S_latente
S_latente <- function(w, mu, sigma) {
  1 - pnorm((log(w) - mu) / sigma)
}


## Proporção de Suscetíveis Observados (Baseado na Série da Sobrevivência)
p_obs_susc <- function(tau, mu, sigma,n) {
    
    # 1. Proteção
    if (tau <= 1e-6 || !is.finite(tau)) return(0)
    
    # 2. Limite Superior da Soma (chão(tau))
    tau_floor <- floor(tau)
    if (tau_floor < 1) return(0) 
    
    # 3. Valores de tempo S(y)
    y_vals <- 0:tau_floor
    S_y <- S_latente(y_vals, mu, sigma)
    
    # 4. P(W > C) (Prob. de Censura)
    P_censura <- (1 / tau_floor) * sum(exp(n * log(pmax(S_y, 1e-100))))
    
    if (!is.finite(P_censura) || is.na(P_censura)) {
         return(0) # Retorna 0 no caso de erro ou NaN
    }
    
    return(P_censura)
}

# Função para encontrar a raiz (tau)
find_tau_root <- function(tau, mu, sigma, p_susc_target) {
    return(p_obs_susc(tau, mu, sigma, SAMPLE_SIZES) - p_susc_target)
}

# =================================================================
# SECTION 2: GERAÇÃO DE DADOS (LND + BN + Censura)
# =================================================================

# Pré-cálculo dos valores de Tau por padrão de covariável
calc_tau_map <- function(params, p_susc_target) {
    cov_patterns <- data.frame(x1 = c(0, 1, 0, 1), x2 = c(0, 0, 1, 1))
    Xe_mu_pat <- as.matrix(cbind(1, cov_patterns$x1, cov_patterns$x2))
    sigma_k <- exp(params[1]) 
    mu_patterns <- Xe_mu_pat %*% params[2:4]
    tau_patterns <- numeric(4)
    
    for (k in 1:4) {
        mu_k <- mu_patterns[k]
        root_result <- tryCatch(
            uniroot(find_tau_root, interval = c(1, 15000), 
                    mu = mu_k, sigma = sigma_k, p_susc_target = p_susc_target,
                    extendInt = "yes",tol = 1e-5)$root,
            error = function(e) NA
        )
        tau_patterns[k] <- root_result
    }
    cov_patterns$tau <- tau_patterns
    return(cov_patterns)
}

# Função principal de geração de dados (Compactada e Vetorizada)
gen_data_lnd_nb <- function(N, params, tau_map) {
    
    x1 <- sample(c(0, 1), size = N, replace = TRUE)
    x2 <- sample(c(0, 1), size = N, replace = TRUE)
    
    # 1. Parâmetros (sigma e phi são escalares)
    Xe <- cbind(1, x1, x2)
    
    sigma <- exp(params[1])
    mu <- Xe %*% params[2:4]
    p0 <- 1 / (1 + exp(-Xe %*% params[5:7]))
    phi <- exp(params[8])
    
    eta <- phi^(-1) * (p0^(-phi) - 1)
    eta[eta < 0] <- 0
    
    # 2. Geração de Riscos (M)
    prob_bn <- 1 / (1 + phi * eta)
    M <- rnbinom(N, size = 1 / phi, prob = prob_bn)
    
    # 3. Geração de Tempos Latentes (W)
    W <- rep(Inf, N)
    idx_susc <- which(M > 0)
    
    for (i in idx_susc) {
        W[i] <- min(rlnd_disc(M[i], mu[i], sigma)) 
    }
    
    # 4. Censura Uniforme Discreta (C)
    tau_i <- tau_map$tau[match(interaction(x1, x2), interaction(tau_map$x1, tau_map$x2))]
    max_tau <- max(tau_map$tau, na.rm = TRUE)
    tau_i[is.na(tau_i) | tau_i <= 0] <- max_tau
    
    C <- ceiling(runif(N, min = 0, max = tau_i))
    
    # 5. Observação (Y) e Status (Delta)
    Y <- pmin(W, C)
    Delta <- as.numeric(Y < Inf & Y == W)
    Y[Y == Inf] <- C[Y == Inf] 
    
    return(data.frame(y = Y, status = Delta, x1 = x1, x2 = x2))
}


# =================================================================
# SECTION 3: FUNÇÕES DE ESTIMAÇÃO ML DISCRETA 
# =================================================================

# Função de Log-Verossimilhança
fvero <- function(theta, data) {
    
    y <- data$y
    status <- data$status
    Xe <- cbind(1, data$x1, data$x2)
    
    # Parâmetros
    sigma <- exp(theta[1]) 
    mu <- Xe %*% theta[2:4]
    p0 <- 1 / (1 + exp(-Xe %*% theta[5:7]))
    phi <- exp(theta[8])
    
    # Proteção e ETA
    phi_safe <- pmax(phi, 1e-6)
    p0_inv_phi <- p0^(-phi_safe) 
    if(any(!is.finite(p0_inv_phi))){ return(-1e20) }
    
    eta <- phi_safe^(-1) * (p0_inv_phi - 1)
    if(any(!is.finite(eta)) || any(eta < 0)){ return(-1e20) }

    # S_LND_fit inline (Reciclagem de sigma)
    S_LND <- function(w, mu_vec) {
        1 - pnorm((log(pmax(w, 1e-6)) - mu_vec) / sigma) 
    }
    
    # Cálculo da Verossimilhança
    vF1 <- 1 - S_LND(y, mu)  # 1 - S_latente(y)
    Spop1 <- (1 + phi * eta * vF1)^(-1 / phi) # S_pop(y)
    
    vF2 <- 1 - S_LND(y + 1, mu)
    Spop2 <- (1 + phi * eta * vF2)^(-1 / phi) # S_pop(y+1)
    
    fpop <- pmax(Spop1 - Spop2, 1e-16) # f_pop(y)
    Spop1 <- pmax(Spop1, 1e-16)
    
    loglik <- sum(status * log(fpop) + (1 - status) * log(Spop1))
    
    if(!is.finite(loglik)){ return(-1e20) }
    
    return(loglik)
}


# otimização
ml_optim <- function(theta_init, data) {
    
   mgg <- optim(theta_init, fvero, data = data, method = "BFGS",
                control = list(fnscale = -1, maxit = 1000))
    
#   theta0 <- theta_init
#   mgg1 <- optim(theta0, fvero, data = data, method = "SANN",
#              control = list(fnscale = -1, maxit = 10000))
#   theta0 <- mgg1$par
#   mgg <- optim(theta0, fvero, data = data, method = "BFGS",
#             control = list(fnscale = -1, maxit = 1000))
    
    if (mgg$convergence != 0) {
        return(list(par = rep(NA, N_PARAM), converged = FALSE))
    }
    
    # Cálculo da Matriz de Informação de Fisher (Hessiana)
    obsinf <- tryCatch(-fdHess(mgg$par, fvero, data = data)$Hessian, error = function(e) NA)
    
    if(any(is.na(obsinf))){
        return(list(par = rep(NA, N_PARAM), converged = FALSE))
    }
    
    covmat <- tryCatch(
        ginv(obsinf),
        error = function(e) matrix(NA, N_PARAM, N_PARAM)
    )
    
    setheta <- sqrt(pmax(0, diag(covmat)))
    
    # Validação rigorosa do erro padrão
    if(any(is.na(setheta)) || any(setheta > 10) || any(setheta <= 0)){
        return(list(par = rep(NA, N_PARAM), converged = FALSE))
    }
    
    return(list(par = mgg$par, se = setheta, converged = TRUE))
}


# =================================================================
# SECTION 4: SIMULAÇÃO MONTE CARLO (EXECUÇÃO E MÉTRICAS)
# =================================================================

TAU_MAP_GLOBAL_DISC <- calc_tau_map(THETA_TRUE, P_SUSC_TARGET)

if(any(is.na(TAU_MAP_GLOBAL_DISC$tau)) || any(TAU_MAP_GLOBAL_DISC$tau <= 0)){
    warning("Erro no cálculo de TAU. Usando fallback.")
    max_tau <- max(TAU_MAP_GLOBAL_DISC$tau, na.rm = TRUE)
    TAU_MAP_GLOBAL_DISC$tau[is.na(TAU_MAP_GLOBAL_DISC$tau) | TAU_MAP_GLOBAL_DISC$tau <= 0] <- max_tau
}

resultados <- list()

for (N in SAMPLE_SIZES) {
    
    cat(sprintf("\n--- Executando simulação DISCRETA (LND) para N = %d (%d repetições) ---\n", N, N_MONTE_CARLO))
    
    estimativas_rep <- matrix(NA, nrow = N_MONTE_CARLO, ncol = N_PARAM, dimnames = list(NULL, names(THETA_TRUE)))
    se_rep <- matrix(NA, nrow = N_MONTE_CARLO, ncol = N_PARAM, dimnames = list(NULL, paste0("SE_", names(THETA_TRUE))))
    
    valid_runs <- 0
    
    for (run in 1:N_MONTE_CARLO) {
        
        set.seed(N + run)
        
        sim_data <- gen_data_lnd_nb(N, THETA_TRUE, TAU_MAP_GLOBAL_DISC)
        ml_result <- ml_optim(THETA_TRUE, sim_data)
        
        if (ml_result$converged) {
            valid_runs <- valid_runs + 1
            estimativas_rep[valid_runs, ] <- ml_result$par
            se_rep[valid_runs, ] <- ml_result$se
        }
        
        if (run %% 50 == 0) cat(sprintf("  Repetição %d concluída. Rodadas válidas: %d\n", run, valid_runs))
    }
    
    # Limpeza e Análise
    estimativas_clean <- estimativas_rep[1:valid_runs, ]
    se_clean <- se_rep[1:valid_runs, ]
    M_valid <- valid_runs
    
    if (M_valid < 10) { warning(sprintf("N=%d: Apenas %d rodadas válidas. Análise limitada.", N, M_valid)) }
    
    # 1. Métrica de Vício e Erro
    bias <- colMeans(estimativas_clean) - THETA_TRUE
    sq_diff <- sweep(estimativas_clean, 2, THETA_TRUE, "-")^2
    rmse <- sqrt(colMeans(sq_diff))
    
    # 2. Métrica de Erro Padrão (SE)
    avg_se <- colMeans(se_clean)
    
    # 3. Métricas do Intervalo de Confiança (IC de 95%)
    Z_95 <- qnorm(0.975)
    
    # AW (Average Width - Largura Média)
    aw <- 2 * Z_95 * avg_se
    
    # CP (Empirical Coverage Probability - Probabilidade de Cobertura Empírica)
    lower_ci <- estimativas_clean - Z_95 * se_clean
    upper_ci <- estimativas_clean + Z_95 * se_clean
    cp <- colMeans(
        sweep(lower_ci, 2, THETA_TRUE, "<") & sweep(upper_ci, 2, THETA_TRUE, ">")
    )
    
    # Armazenamento dos resultados
    resultados[[as.character(N)]] <- data.frame(
        N = N,
        Parameter = names(THETA_TRUE),
        True_Value = THETA_TRUE,
        Bias = bias,
        RMSE = rmse,
        Avg_SE = avg_se,
        AW = aw,      
        CP = cp,      
        Valid_Runs = M_valid
    )
    cat(sprintf("- N=%d: %d de %d rodadas válidas (Taxa de Convergência: %.1f%%).\n", 
                N, M_valid, N_MONTE_CARLO, (M_valid / N_MONTE_CARLO) * 100))
}

# =================================================================
# SECTION 5: ANÁLISE GRÁFICA
# =================================================================

cat("\n--- GERANDO GRÁFICOS DE DESEMPENHO ---\n")

df_resultados <- do.call(rbind, resultados)

# Reestrutura os dados para ggplot (long format)
df_long <- pivot_longer(df_resultados, 
                        cols = c(Bias, RMSE, Avg_SE, AW, CP),
                        names_to = "Metric",
                        values_to = "Value")

# 5.1. Gráfico de Vício (Bias)
plot_bias <- ggplot(subset(df_long, Metric == "Bias"), 
                    aes(x = N, y = Value, color = Parameter)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~ Parameter, scales = "free_y") +
    labs(title = "Vício (Bias) vs. Tamanho Amostral (N)",
         y = "Vício (Bias)", x = "Tamanho Amostral (N)") +
    theme_minimal() +
    theme(legend.position = "none")

# 5.2. Gráfico de Precisão e Largura (RMSE, Avg_SE e AW)
df_precision <- subset(df_long, Metric %in% c("RMSE", "Avg_SE", "AW"))

plot_precision <- ggplot(df_precision, aes(x = N, y = Value, color = Metric, linetype = Metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ Parameter, scales = "free_y") +
    labs(title = "Precisão e Largura Média (RMSE, Avg. SE, AW) vs. N",
         y = "Valor da Métrica", x = "Tamanho Amostral (N)", color = "Métrica", linetype = "Métrica") +
    scale_color_manual(values = c("Avg_SE" = "blue", "AW" = "darkgreen", "RMSE" = "red")) +
    scale_linetype_manual(values = c("Avg_SE" = "solid", "AW" = "dashed", "RMSE" = "dotted")) +
    theme_minimal()

# 5.3. Gráfico de Probabilidade de Cobertura (CP)
plot_coverage <- ggplot(subset(df_long, Metric == "CP"), 
                       aes(x = N, y = Value, color = Parameter)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    facet_wrap(~ Parameter) +
    labs(title = "Probabilidade de Cobertura Empírica (CP) vs. Tamanho Amostral (N)",
         y = "Probabilidade de Cobertura", x = "Tamanho Amostral (N)") +
    ylim(0.8, 1.0) + 
    theme_minimal() +
    theme(legend.position = "none")

# Exibição dos gráficos
print(plot_bias)
print(plot_precision)
print(plot_coverage)

# =================================================================
# SECTION 6: EXIBIÇÃO DOS RESULTADOS TABULADOS
# =================================================================

print("--- RESULTADOS DA SIMULAÇÃO MONTE CARLO (LND Discreta - CÓDIGO REFATORADO) ---")
# Exibe a tabela completa como um dataframe
print(df_resultados)
