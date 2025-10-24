# Install and load packages
install.packages("FrF2")
library(FrF2)

# Type 1 error experiment

# Factor A = Sample size
fA <- function(levelA, n) {
  if (levelA == -1) { # Low level
    data <- rnorm(n = 5, mean = 0, sd = 1)
  } else if (levelA == 1) { # High level
    data <- rnorm(n = 300, mean = 0, sd = 1)
  }
  return(data)
}

# Factor B = Distribution
fB <- function(levelB, data) {
  if (levelB == -1) { # Normal distribution
    transformed_data <- data
  } else if (levelB == 1) { # t-distribution
    transformed_data <- rt(length(data), df = 3)
  }
  return(transformed_data)
}

# Factor C = Test type
fC <- function(levelC, data) {
  p_values <- c()
  for (i in 1:100) {
    sample1 <- sample(data, length(data) / 2)
    sample2 <- sample(data, length(data) / 2)

    if (levelC == -1) { # t-test
      test <- t.test(sample1, sample2)
    } else if (levelC == 1) { # z-test (assuming known variance)
      term_sqrt <- var(sample1) / length(sample1) + var(sample2) / length(sample2)
      test_statistic <- (mean(sample1) - mean(sample2)) / sqrt(term_sqrt)
      p_value <- 2 * (1 - pnorm(abs(test_statistic)))
      test <- list(p.value = p_value)
    }

    p_values[i] <- test$p.value
  }
  error_rate <- sum(p_values < 0.05) / 100
  return(error_rate)
}

# 2^3 factorial experiment
experiment <- function(n, r) {
  data <- fA(r[1], n)
  transformed_data <- fB(r[2], data)
  error_rate <- fC(r[3], transformed_data)
  return(error_rate)
}

# Design matrix
X <- cbind(
  c(-1, -1, -1, -1, +1, +1, +1, +1),  # A
  c(-1, -1, +1, +1, -1, -1, +1, +1),  # B
  c(-1, +1, -1, +1, -1, +1, -1, +1)   # C
)
colnames(X) <- c("A", "B", "C")

num_replicates <- 10
X_full <- X
for (i in 2:num_replicates) {
  X_full <- rbind(X_full, X)
}
colnames(X_full) <- c("A", "B", "C")

# One run for plots
random_order <- sample(1:nrow(X_full))
result <- c()

for (i in random_order) {
  data <- fA(X_full[i, "A"])
  transformed_data <- fB(X_full[i, "B"], data)
  type1error <- fC(X_full[i, "C"], transformed_data)
  result <- c(result, type1error)
}

dataset <- data.frame(X_full[random_order, ], Response = result)

# Fit model
model <- lm(Response ~ A * B * C, data = dataset)
summary(model)
DanielPlot(model)
MEPlot(model)
IAPlot(model)

# --- SIMULATION LOOP ---
num_runs <- 100
coef_list <- list()
pval_list <- list()
model_pvals <- numeric(num_runs)
r_sq_vec <- numeric(num_runs)
adj_r_sq_vec <- numeric(num_runs)
resid_se_vec <- numeric(num_runs)
residual_list <- list()

for (i in 1:num_runs) {
  result <- c()
  random_order <- sample(1:nrow(X_full))
  for (j in random_order) {
    data <- fA(X_full[j, "A"])
    transformed_data <- fB(X_full[j, "B"], data)
    type1error <- fC(X_full[j, "C"], transformed_data)
    result <- c(result, type1error)
  }

  dataset <- data.frame(X_full[random_order, ], Response = result)
  model <- lm(Response ~ A * B * C, data = dataset)
  summary_model <- summary(model)

  coef_list[[i]] <- coef(summary_model)[, "Estimate"]
  pval_list[[i]] <- coef(summary_model)[, "Pr(>|t|)"]
  r_sq_vec[i] <- summary_model$r.squared
  adj_r_sq_vec[i] <- summary_model$adj.r.squared
  resid_se_vec[i] <- summary_model$sigma
  residual_list[[i]] <- residuals(model)

  f_stat <- summary_model$fstatistic
  p_model <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
  model_pvals[i] <- p_model
}

# --- ANALYSIS ---
coef_df <- do.call(rbind, coef_list)
pval_df <- do.call(rbind, pval_list)

coef_stats <- apply(coef_df, 2, function(x) c(mean = mean(x), sd = sd(x)))
sig_props <- apply(pval_df, 2, function(p) mean(p < 0.05))

summary_table <- data.frame(
  Term = colnames(coef_df),
  Mean_Coefficient = coef_stats["mean", ],
  SD_Coefficient = coef_stats["sd", ],
  Prop_Significant_p_lt_0.05 = sig_props
)

print(summary_table)

# --- PLOTS ---
barplot(
  summary_table$Prop_Significant_p_lt_0.05,
  names.arg = summary_table$Term,
  las = 2,
  col = "skyblue",
  ylim = c(0, 1),
  main = "Proportion of p-values < 0.05 per Term",
  ylab = "Proportion Significant"
)

boxplot(
  coef_df,
  main = "Distribution of Coefficient Estimates",
  ylab = "Estimate",
  las = 2,
  col = "lightgreen"
)
abline(h = 0, col = "gray")

hist(
  model_pvals,
  main = "Distribution of Model P-values",
  xlab = "Model p-value",
  col = "lightblue",
  breaks = 20
)
abline(v = 0.05, col = "red", lty = 2)

summary_table_no_intercept <- subset(summary_table, Term != "(Intercept)")
barplot(
  summary_table_no_intercept$Mean_Coefficient,
  names.arg = summary_table_no_intercept$Term,
  las = 2,
  col = "steelblue",
  main = "Average Effect Size Across Simulations",
  ylab = "Mean Coefficient Estimate"
)
abline(h = 0, col = "gray", lty = 2)

# --- QQ PLOTS ---
par(mfrow = c(3, 3))
for (i in 1:9) {
  qqnorm(residual_list[[i]], main = paste("Q-Q Plot Run", i))
  qqline(residual_list[[i]], col = "red")
}
par(mfrow = c(1, 1))
