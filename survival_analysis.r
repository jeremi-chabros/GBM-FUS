library(survival)
library(lubridate)
library(ggplot2)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(rms)

method <- "CEM"

if (method == "optimal") {
    df <- read.csv("matched_data_optimal.csv")
} else {
    df <- read.csv("matched_data.csv")
    method <- "CEM"
}

# Convert dates to proper format
df$SurgeryDate <- dmy(df$SurgeryDate)
df$DeathCensorDate <- dmy(df$DeathCensorDate)

# Survival defined in number of months since surgery
df$Survival <- as.duration(df$SurgeryDate %--% df$DeathCensorDate) / dmonths(1)

# Change names for plotting
df$FUS <- ifelse(df$FUS == 1, "FUS", "Standard of care")

summary(df)


# Fit survival curve
fit_surv <- survfit2(Surv(df$Survival, df$Dead) ~ FUS, data = df)
ggsurv_plot <- fit_surv %>%
    ggsurvfit() +
    labs(
        x = "Time [Months]",
        y = "Overall survival probability"
    ) +
    add_confidence_interval() +
    add_risktable()

# Fitting the survival curve
fit_surv <- survfit(Surv(df$Survival, df$Dead) ~ FUS, data = df)
fit_surv

# Save the plot
ggsave(paste("plots/survival_curve_", method, ".png", sep=""), ggsurv_plot, width = 10, height = 7, res=300)

# Only 33 observations with rule of 10 gives 3 covariates to avoid overfitting
# res.cox = coxph(Surv(df$Survival, df$Dead) ~ FUS + Age + Gender + Race + MGMT + IDH + X1p19q, data = df)
res.cox <- coxph(Surv(df$Survival, df$Dead) ~ FUS + IDH, data = df)
summary(res.cox)
res.cox

cox.zph_fit <- cox.zph(res.cox)
print(cox.zph_fit)

# Plotting
png(paste("plots/Schoenfeld_residuals_", method, ".png", sep=""), units="in", width=5, height=5, res=300)
plot(cox.zph_fit)
dev.off()



### Fit Accelerated Failure Time (AFT) Model
library(flexsurv)

fit_weibull <- flexsurvreg(Surv(Survival, Dead) ~ FUS, data = df, dist = "weibull")
fit_exp <- flexsurvreg(Surv(Survival, Dead) ~ FUS, data = df, dist = "exp")
fit_gamma <- flexsurvreg(Surv(Survival, Dead) ~ FUS, data = df, dist = "gamma")

png(paste("plots/residuals_weibull_qqplot_", method, ".png", sep=""), units="in", width=5, height=5, res=300)
residuals_weibull <- residuals(fit_weibull, type = "response")
qqnorm(residuals_weibull)
qqline(residuals_weibull)
dev.off()

png(paste("plots/residuals_exp_qqplot_", method, ".png", sep=""), units="in", width=5, height=5, res=300)
residuals_exp <- residuals(fit_exp, type = "response")
qqnorm(residuals_exp)
qqline(residuals_exp)
dev.off()

png(paste("plots/residuals_gamma_qqplot_", method, ".png", sep=""), units="in", width=5, height=5, res=300)
residuals_gamma <- residuals(fit_gamma, type = "response")
qqnorm(residuals_gamma)
qqline(residuals_gamma)
dev.off()

# Compare
aic_bic_data <- data.frame(
    Distribution = c("Weibull", "Exponential", "Gamma"),
    AIC = c(AIC(fit_weibull), AIC(fit_exp), AIC(fit_gamma)),
    BIC = c(BIC(fit_weibull), BIC(fit_exp), BIC(fit_gamma))
)
print(aic_bic_data)

fit_aft <- flexsurvreg(Surv(Survival, Dead) ~ FUS + MGMT + IDH, data = df, dist = "weibull")
fit_aft

# Function for Sample Size Calculation for Log-rank Test
calculate_sample_size <- function(alpha, power, hazard_ratio) {
    # Z-scores for alpha and power
    z_alpha <- qnorm(1 - alpha / 2)
    z_power <- qnorm(1 - power)

    # Calculate sample size for each group
    n_each_group <- ((z_alpha + z_power)^2) / (log(hazard_ratio)^2)

    return(n_each_group)
}

# Parameters
alpha <- 0.05 # Significance level
power <- 0.8 # Power
hazard_ratio <- 1.76218 # Hazard ratio between the two groups

# Sample Size Calculation
n <- calculate_sample_size(alpha, power, hazard_ratio)

# Output the calculated sample size for each group
print(paste("Sample size needed for each group: ", round(n)))
