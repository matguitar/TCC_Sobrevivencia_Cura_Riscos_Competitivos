library(survival)
library(nlme)    # Para fdHess (Hessiana)
library(MASS)    # Para ginv (Pseudo-inversa)
library(stats)    # Para uniroot, pnorm, qnorm, optim, rnbinom
library(ggplot2) # Para gráficos
library(tidyr)    # Para pivot_longer (processamento de dados para gráficos)

set.seed(1)

# =================================================================
# SECTION 0: DEFINIÇÕES E CONFIGURAÇÃO (AJUSTADO PARA 6 PARÂMETROS)
# =================================================================

logit <- function(p) {
    log(p / (1 - p))
}

# Parâmetros verdadeiro.5s (True Values) 
PARAM_FIXOS_NOVO <- list(
    # 1. Parâmetro de Dispersão LND
    beta_sigma_int = 0.2,         # log(sigma)
    
    # 2. Parâmetro de Mu (Lognormal) - APENAS INTERCEPTO
    beta_mu_int = 5,          # Intercepto (Nível z=1, z=2, z=3)
    
    # 3-5. Parâmetros de P0 (Prob. de cura) 
    beta_p0_int = logit(0.20), # Intercepto (Nível z=1)
    beta_p0_z2 = logit(0.35), 
    beta_p0_z3 = logit(0.50), 
    
    # 6. Parâmetro de Dispersão NB
    beta_phi_int = 0.8        # log(phi)
)

THETA_TRUE <- unlist(PARAM_FIXOS_NOVO)
N_PARAM <- length(THETA_TRUE) # Deve ser 6
names(THETA_TRUE) <- names(PARAM_FIXOS_NOVO)

# Variáveis do estudo de simulação
SAMPLE_SIZES <- c(100, 200, 300, 400, 500, 700, 1000)
N_MONTE_CARLO <- 1000
P_SUSC_TARGET <- 0.05
P_OBS_TARGET <- 1 - P_SUSC_TARGET 


# =================================================================
# SECTION 1: FUNÇÕES LND DISCRETA 
# =================================================================

# Função de Sobrevivência Lognormal Discreta
sf_lnd_disc <- function(w, mu, sigma) {
    w <- pmax(w, 1e-6) 
    1 - pnorm((log(w) - mu) / sigma)
}

# Geração de Tempo Aleatório LND Discreto
rlnd_disc <- function(n, mu, sigma) {
    U <- runif(n)
    W <- ceiling(exp(mu + sigma * qnorm(1 - U))) - 1
    pmax(0, W)
}

## Proporção de Suscetíveis Observados
p_obs_susc <- function(tau, mu, sigma) {
    
    if (tau <= 1e-6 || !is.finite(tau)) return(0)
    
    tau_floor <- floor(tau)
    if (tau_floor < 1) return(0) 
    
    y_vals <- 1:tau_floor
    S_y <- sf_lnd_disc(y_vals, mu, sigma)
    
    P_censura <- (1 / tau) * sum(S_y)
    
    P_observada <- pmax(0, 1 - P_censura)
    
    if (!is.finite(P_observada) || is.na(P_observada)) {
          return(0) 
    }
    
    return(P_observada)
}

# Função para encontrar a raiz (tau)
find_tau_root <- function(tau, mu, sigma, p_obs_target) {
    return(p_obs_susc(tau, mu, sigma) - p_obs_target)
}

# =================================================================
# SECTION 2: GERAÇÃO DE DADOS (AJUSTADO PARA MU CONSTANTE)
# =================================================================

## Pré-cálculo dos valores de Tau por padrão de covariável (APENAS MU)
calc_tau_map <- function(params, p_obs_target) {
    # 1. Definir a única covariável de 3 níveis
    cov_patterns <- data.frame(z = c(1, 2, 3))
    
    # 2. Criar variáveis dummy (necessário apenas para a geração de dados, mas aqui não importa para mu)
    cov_patterns$dummy_2 <- ifelse(cov_patterns$z == 2, 1, 0)
    cov_patterns$dummy_3 <- ifelse(cov_patterns$z == 3, 1, 0)

    # 3. Matriz de Desenho (APENAS INTERCEPTO PARA MU)
    # params[1] = log(sigma_k), params[2] = beta_mu_int
    Xe_mu_pat <- as.matrix(1) # Matriz com apenas 1 (Intercepto)
    
    sigma_k <- exp(params[1]) 
    
    # Cálculo dos valores de mu: mu é constante (usa apenas params[2])
    mu_constant <- params[2]
    
    # 4. Cálculo do Tau: Tau é constante, pois mu e sigma são constantes
    tau_patterns <- numeric(3) 
    
    # O Tau será o mesmo para os 3 níveis
    root_result <- tryCatch(
        uniroot(find_tau_root, interval = c(1e-6, 500), 
                        mu = mu_constant, sigma = sigma_k, p_obs_target = p_obs_target,
                        extendInt = "yes")$root,
        error = function(e) NA
    )
    tau_patterns[] <- root_result # Aplica o mesmo Tau a todos os 3 níveis
    
    cov_patterns$tau <- tau_patterns
    return(cov_patterns)
}

## Função principal de geração de dados (NOVO CENÁRIO)
gen_data_lnd_nb <- function(N, params, tau_map) {
    
    # 1. Geração da Covariável Única (z) com 3 níveis
    z <- sample(c(1, 2, 3), size = N, replace = TRUE)
    
    # Criação das variáveis dummy (para P0)
    dummy_2 <- ifelse(z == 2, 1, 0)
    dummy_3 <- ifelse(z == 3, 1, 0)
    
    # 2. Parâmetros do Modelo (sigma, mu, p0, phi)
    sigma <- exp(params[1])
    
    # MU É CONSTANTE (usa apenas params[2])
    mu <- rep(params[2], N)
    
    # P0 USA AS COVARIÁVEIS (params 3, 4, 5)
    Xe_p0 <- cbind(1, dummy_2, dummy_3)
    p0_lin_pred <- Xe_p0 %*% params[3:5] # Posições 3, 4 e 5 no vetor params
    p0 <- 1 / (1 + exp(-p0_lin_pred))
    
    # PHI É CONSTANTE (params[6])
    phi <- exp(params[6])
    
    # Parâmetro eta para a Binomial Negativa
    eta <- phi^(-1) * (p0^(-phi) - 1)
    eta[eta < 0] <- 0
    
    # 3. Geração de Riscos (M) usando a Binomial Negativa
    prob_bn <- 1 / (1 + phi * eta)
    M <- rnbinom(N, size = 1 / phi, prob = prob_bn)
    
    # 4. Geração de Tempos Latentes (W)
    W <- rep(Inf, N)
    idx_susc <- which(M > 0)
    
    # O rlnd_disc é vetorizado, mas o mu[i] tem que ser passado como valor único
    for (i in idx_susc) {
        # Como mu é constante, poderíamos usar: W[idx_susc] <- rlnd_disc(M[idx_susc], mu[1], sigma)
        # Mas para manter a lógica de min(rlnd_disc(M[i])), mantemos o loop:
        W[i] <- min(rlnd_disc(M[i], mu[i], sigma)) 
    }
    
    # 5. Censura Uniforme Discreta (C)
    tau_i <- tau_map$tau[match(z, tau_map$z)] # Tau ainda é constante, mas esta linha o recupera
    
    max_tau <- max(tau_map$tau, na.rm = TRUE)
    tau_i[is.na(tau_i) | tau_i <= 0] <- max_tau
    
    C <- ceiling(runif(N, min = 0, max = tau_i))
    
    # 6. Observação (Y) e Status (Delta)
    Y <- pmin(W, C)
    Delta <- as.numeric(Y < Inf & Y == W) 
    Y[Y == Inf] <- C[Y == Inf] 
    
    return(data.frame(y = Y, status = Delta, z = z))
}

# =================================================================
# SECTION 3: FUNÇÕES DE ESTIMAÇÃO ML DISCRETA (AJUSTADO PARA 6 PARÂMETROS)
# =================================================================


fvero <- function(theta, data) {
    
    y <- data$y
    status <- data$status 
    
    # 1. Criação das variáveis dummy a partir da covariável 'z'
    z <- data$z
    dummy_2 <- ifelse(z == 2, 1, 0)
    dummy_3 <- ifelse(z == 3, 1, 0)
    
    # 2. Matriz de Desenho para P0 (Intercepto + dummy_2 + dummy_3)
    Xe_p0 <- cbind(1, dummy_2, dummy_3)
    
    # Theta deve ter 6 parâmetros: log(sigma), mu, 3 para p0, log(phi)
    
    # Parâmetros
    sigma <- exp(theta[1]) 
    
    # MU É CONSTANTE (usa apenas theta[2])
    mu <- rep(theta[2], length(y))
    
    # P0 USA AS COVARIÁVEIS (theta[3:5])
    p0_lin_pred <- Xe_p0 %*% theta[3:5] 
    p0 <- 1 / (1 + exp(-p0_lin_pred))
    
    # PHI É CONSTANTE (theta[6])
    phi <- exp(theta[6])
    
    # Proteção e ETA
    phi_safe <- pmax(phi, 1e-6)
    p0_inv_phi <- p0^(-phi_safe) 
    if(any(!is.finite(p0_inv_phi))){ return(-1e20) }
    
    eta <- phi_safe^(-1) * (p0_inv_phi - 1)
    if(any(!is.finite(eta)) || any(eta < 0)){ return(-1e20) }

    # S_LND_fit inline (Reciclagem de sigma)
    # S_LND é a função de Sobrevivência (Survival) da distribuição Lognormal Discreta (S_latente)
    S_LND <- function(w, mu_vec) {
        # Usamos o mu constante (theta[2])
        1 - pnorm((log(pmax(w, 1e-6)) - mu_vec) / sigma) 
    }
    
    # Cálculo da Verossimilhança do modelo LND-NB
    
    # S_pop(y) (Probabilidade de Sobreviver a y - Inclusiva)
    vF1 <- 1 - S_LND(y, mu)  # 1 - S_latente(y)
    Spop1 <- (1 + phi * eta * vF1)^(-1 / phi) 
    
    # S_pop(y+1) (Probabilidade de Sobreviver a y+1 - Exclusiva)
    vF2 <- 1 - S_LND(y + 1, mu)
    Spop2 <- (1 + phi * eta * vF2)^(-1 / phi) 
    
    # f_pop(y) = P(Y=y) = S_pop(y) - S_pop(y+1)
    fpop <- pmax(Spop1 - Spop2, 1e-16) 
    Spop1 <- pmax(Spop1, 1e-16)
    
    # Log-Verossimilhança:
    # Se status=1 (evento), contribui com log(fpop)
    # Se status=0 (censura), contribui com log(Spop(y))
    loglik <- sum(status * log(fpop) + (1 - status) * log(Spop1))
    
    if(!is.finite(loglik)){ return(-1e20) }
    
    return(loglik)
}

# otimização
ml_optim <- function(theta_init, data) {
    

   mgg <- optim(theta_init, fvero, data = data, method = "BFGS",
             control = list(fnscale = -1, maxit = 1000))
    
    
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

# Recalcula o TAU_MAP com os NOVOS parâmetros (Apenas 1 para mu)
TAU_MAP_GLOBAL_DISC <- calc_tau_map(THETA_TRUE, P_OBS_TARGET)

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
        AW = aw,
        CP = cp     
        
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
                        cols = c(Bias, RMSE, AW, CP), 
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

# 5.2. Gráfico de Precisão e Largura (RMSE, e AW)
df_precision <- subset(df_long, Metric %in% c("RMSE", "AW"))

plot_precision <- ggplot(df_precision, aes(x = N, y = Value, color = Metric, linetype = Metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ Parameter, scales = "free_y") +
    labs(title = "Precisão e Largura Média (RMSE, Avg. SE, AW) vs. N",
         y = "Valor da Métrica", x = "Tamanho Amostral (N)", color = "Métrica", linetype = "Métrica") +
    scale_color_manual(values = c("AW" = "darkgreen", "RMSE" = "blue")) +
    scale_linetype_manual(values = c("AW" = "dashed", "RMSE" = "dotted")) +
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

print("--- RESULTADOS DA SIMULAÇÃO MONTE CARLO (LND Discreta - 6 PARÂMETROS REFATORADO) ---")
# Exibe a tabela completa como um dataframe
print(df_resultados)
