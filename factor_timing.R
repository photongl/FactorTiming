library(nloptr)
library(readr)
library(lubridate)
library(PerformanceAnalytics)


# ======== Parameters ========
returns_data <- "return_data.csv"
factor_data <- "factor_data.csv"
signal_data <- "signal_data_momentum.csv"
transaction_cost_data <- "transaction_cost_data.csv"
ref_portfolio_data <- "ref_portfolio_data.csv"

REF_PF_WEIGHTS <- c(0.6, 0.4)
MAX_ACTIVE_RISK <- 0.15   # annualized
BUDGET <- 1.0         # pf weights sum up to this value
IC <- 0.5             # information coefficient: how correlated our signal is with realized returns (essentially confidence in signal)
lambda <- 2           # risk-aversion coefficient
theta <- 0            # transaction-cost coefficient, set to zero for no costs
vcov_window = 131     # (in days) - window used to calculate factor vcov matrix
beta_window = 131     # (in days) - window used to calculate asset-beta against factors (ie factor loadings)
periods_in_year = 262 # set this to 12, if data-frequency is monthly, 262 if daily etc.
max_optimizer_iterations = 10000
beta_significance_threshold = 0.995 # (1 - two-tailed type-I error probability) for factor-loading significance
MAX_SHORT_SIZE = 0.0

# ========= Helper functions ==========
utility <- function(h, alpha, beta, h_ref, beta_ref, V, lambda, theta, tcost, h_i)
{
  u = (-1.0 * ((t(h) %*% t(beta) - t(h_ref) %*% t(beta_ref)) %*% alpha 
               - lambda * ((t(h) %*% t(beta) - t(h_ref) %*% t(beta_ref)) %*% V %*% (beta %*% h - beta_ref %*% h_ref)) 
               - theta * (tcost %*% abs(h - h_i))))
  return(u)
}

utility_grad <- function(h, alpha, beta, h_ref, beta_ref, V, lambda, theta, tcost, h_i)
{
  -(t(beta) %*% alpha 
    - 2.0 * lambda * (t(beta) %*% V %*% (beta %*% h - beta_ref %*% h_ref)) 
    - theta * (t(tcost) * sign(h - h_i)))
}

budget_constraint <- function(h, alpha, beta, h_ref, beta_ref, V, lambda, theta, tcost, h_i)
{
  return(sum(h) - BUDGET)
}

grad_budget_constraint <- function(h, alpha, beta, h_ref, beta_ref, V, lambda, theta, tcost, h_i)
{
  return(rep(1, num_assets))
}

risk_constraint <- function(h, alpha, beta, h_ref, beta_ref, V, lambda, theta, tcost, h_i)
{
  return(((t(h) %*% t(beta) - t(h_ref) %*% t(beta_ref)) %*% V %*% (beta %*% h - beta_ref%*%h_ref))  - MAX_ACTIVE_RISK * MAX_ACTIVE_RISK)
}

grad_risk_constraint <- function(h, alpha, beta, h_ref, beta_ref, V, lambda, theta, tcost, h_i)
{
  return(2 * (t(beta) %*% V %*% (beta %*% h - beta_ref %*% h_ref)))
}

beta_matrix <- function(asset_returns, factors, rf, t_critical = 2.594724)
{
  beta <- matrix(0, ncol = ncol(asset_returns), nrow = ncol(factors) + 1)
  for(n in 1:ncol(asset_returns))
  {
    y <- asset_returns[ , n] - rf
    fit <- lm(y ~ ., data = factors)
    significant <- ifelse(abs(summary(fit)$coefficients[ , 3]) > t_critical, 1, 0)
    significant[is.na(significant)] <- 0
    beta[ , n] <- fit$coef * significant
  }
  rownames(beta) <- c("alpha", colnames(factors))
  colnames(beta) <- c(colnames(asset_returns))
  return(beta)
}

end_of_month_flag <- function(date_vector)
{
  month_date_vector <- month(date_vector)
  month_date_vector_1 <- c(month_date_vector[2:length(month_date_vector)]
                           , month_date_vector[length(month_date_vector)])
  
  return((month_date_vector_1 - month_date_vector)!=0)
}

# ======== Main Analysis ========
# read data
df <- read_csv(returns_data)
df_factors <- read_csv(factor_data)
df_signals <- read_csv(signal_data)
df_transaction_costs <- read_csv(transaction_cost_data)
df_ref_portfolio <- read_csv(ref_portfolio_data)
rf <- df$`H15T1M Index`

# data size and attributes
num_assets <- ncol(df) - 1
num_factors <- ncol(df_factors) - 1
 
asset_labels <- colnames(df[, -1])
factor_labels <- colnames(df_factors[, -1])
tcost <- t(as.matrix(df_transaction_costs[, -1]))
t_stat_threshold = qt(beta_significance_threshold, beta_window - num_factors)

date_vector <- as.Date(df$Date)
eom_vector <- end_of_month_flag(date_vector)

N <- length(date_vector)
h_matrix <- matrix(nrow = N, ncol = num_assets)
matrix_beta <- matrix(nrow = length(eom_vector[eom_vector==TRUE]), ncol = num_factors + 1)
matrix_beta_ref <- matrix(nrow = length(eom_vector[eom_vector==TRUE]), ncol = num_factors + 1)
colnames(matrix_beta) <- c("alpha", factor_labels)
colnames(matrix_beta_ref) <- c("alpha", factor_labels)
colnames(h_matrix) <- asset_labels
pf_return <- numeric(N)
gross_pf_return <- numeric(N)

# optimizer options
set.seed(5)
h_i <- rep(BUDGET / num_assets, num_assets)

lb <- rep(-MAX_SHORT_SIZE, num_assets)
ub <- rep(100, num_assets)
local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel" = 1.0e-16 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG_EQ",
              "xtol_rel" = 1.0e-16,
              "maxeval" = max_optimizer_iterations,
              "local_opts" = local_opts )

# reference portfolio returns
ref_pf_r <- as.matrix(df_ref_portfolio[, -1]) %*% REF_PF_WEIGHTS

# start optimization
window_size <- max(vcov_window, beta_window)
r_factor_matrix <- as.matrix(cbind(const = 1, df_factors[, -1]))
k = 1
for (i in (window_size + 1):N)
{
  # get asset returns
  r_asset_period <- as.matrix(df[(i+1), -1])
  
  h <- h_i
  # optimize portfolio at month-end
  if (eom_vector[i]==TRUE)
  {
    # get factor returns
    r_factor <- r_factor_matrix[(i - vcov_window):i, ]
    # get covariance matrix
    V <- cov(r_factor) * sqrt(periods_in_year)
    
    # get optimal portfolio
    h_ref <- REF_PF_WEIGHTS
    alpha <- t(as.matrix(cbind(const=1, df_signals[df_signals$Date==date_vector[i], -1]))) * IC
    beta <- beta_matrix(df[((i-beta_window):i), -1]
                        , df_factors[((i-beta_window):i), -1]
                        , rf = rf[((i-beta_window):i)]
                        , t_critical = t_stat_threshold)
    beta_ref <- beta_matrix(df_ref_portfolio[((i-beta_window):i), -1]
                            , df_factors[((i-beta_window):i), -1]
                            , rf = rf[((i-beta_window):i)]
                            , t_critical = t_stat_threshold)
    
    res <- nloptr( x0=h,
                   eval_f=utility,
                   eval_grad_f = utility_grad,
                   lb=lb,
                   ub=ub,
                   eval_g_ineq=risk_constraint,
                   eval_g_eq=budget_constraint,
                   eval_jac_g_ineq = grad_risk_constraint,
                   eval_jac_g_eq = grad_budget_constraint,
                   opts=opts,
                   alpha=alpha, beta=beta, V=V, lambda=lambda, theta=theta, tcost=tcost, h_i=h_i,
                   h_ref = h_ref, beta_ref = beta_ref)
    
    h <- res$solution
    matrix_beta[k, ] <- t(beta %*% h)
    matrix_beta_ref[k, ] <- t(beta_ref %*% h_ref)
    cat(sprintf("Period %s: %s\tSum (pf wts): %s, Active Risk: %s\nOpt status: %s\nIterations: %s\n"
                , i
                , date_vector[i]
                , sum(h)
                , sqrt(((t(h) %*% t(beta) - t(h_ref) %*% t(beta_ref)) %*% V %*% (beta %*% h - beta_ref%*%h_ref)))
                , res$message
                , res$iterations))
    k = k + 1
  }
  
  # store portfolio weights and returns
  h_matrix[i,] <- h
  pf_return[i] <- r_asset_period %*% h - tcost %*% abs(h - h_i)
  gross_pf_return[i] <- r_asset_period %*% h
  h_i <- h
}

# create portfolio-return timeseries
pf_r <- xts(pf_return, order.by = date_vector)
gross_pf_r <- xts(gross_pf_return, order.by = date_vector)
plot_cols <- sample(colors(distinct = TRUE), num_assets)
basket <- xts(x = h_matrix[eom_vector, ], order.by = date_vector[eom_vector])
beta_timeseries <- xts(x = matrix_beta, order.by = date_vector[eom_vector])
beta_ref_timeseries <- xts(x = matrix_beta_ref, order.by = date_vector[eom_vector])
returns <- xts(x = cbind(pf_r, gross_pf_r, ref_pf_r), order.by = as.Date(date_vector))[(window_size):(N-1)]
rf <- xts(x = cbind(rf=rf), order.by = as.Date(date_vector))[(window_size):(N-1)]
colnames(returns) <- c("Active Portfolio (Net)", "Active Portfolio (Gross)", "Reference Portfolio")

# ========== Display results ============
charts.PerformanceSummary(returns, wealth.index = TRUE, main = "Portfolio Performance", colorset=rainbow8equal, cex.axis=1.25, cex.legend=1.5)
layout(matrix(1))
plot.zoo(x=basket, ylab="weight", main = "Portfolio weights", screen=1, col=plot_cols, lwd=2)
legend(x = "topleft", legend = asset_labels, lty = 1, col = plot_cols, lwd=2)
plot.zoo(x=beta_timeseries, ylab="weight", main = "Active Portfolio Beta Evolution", screen=1, col=plot_cols, lwd=2)
legend(x = "topleft", legend = factor_labels, lty = 1, col = plot_cols, lwd=2)
plot.zoo(x=beta_ref_timeseries, ylab="weight", main = "Reference Portfolio Beta Evolution", screen=1, col=plot_cols, lwd=2)
legend(x = "topleft", legend = factor_labels, lty = 1, col = plot_cols, lwd=2)
hist(returns[,1], breaks=100, main="Daily Returns Distribution")

cat("========== Performance metrics ============\n")
print(SharpeRatio.annualized(returns, rf))
print(StdDev.annualized(returns))
print(Return.annualized(returns))
print(skewness(returns))
cat(sprintf("Max drawdown (active-net): %s\n", maxDrawdown(returns[,1])))
cat(sprintf("Max drawdown (active-gross): %s\n", maxDrawdown(returns[,2])))
cat(sprintf("Max drawdown (reference): %s\n", maxDrawdown(returns[,3])))

# write portfolio weights and returns to file
write.zoo(basket, "output_portfolio_weights_momentum.csv", quote=FALSE, sep=",")
write.zoo(returns, "output_pf_returns_momentum.csv", quote=FALSE, sep=",")
write.zoo(beta_timeseries, "output_beta_momentum.csv", quote=FALSE, sep=",")
write.zoo(beta_ref_timeseries, "output_beta_ref_momentum.csv", quote=FALSE, sep=",")
