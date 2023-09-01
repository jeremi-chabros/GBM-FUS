library(MatchIt)
library("marginaleffects")
library(survival)
library(lubridate)

# Read data
read_csv <- read.csv("new.csv")
df <- na.omit(read_csv) # Drop missing values
df <- df[!(is.na(df$DeathCensorDate) | df$DeathCensorDate==""), ]
df <- df[!(df$X1p19q==1), ]
df <- df[df$Chemotherapy==1,]
df <- df[df$Radiotherapy==1,]
df <- df[!(df$TumorSize=="Missing"),]

# Dichotomize race
df$Race <- ifelse(df$Race == "White", "White", "Non-white")
# mymodel <- as.formula(FUS ~ Age + Gender + Race + MGMT + IDH + X1p19q)
mymodel <- as.formula(FUS ~ Age + Gender + Race + TumorSize)

m.out0 <- matchit(mymodel,
    data = df,
    method = NULL, distance = "glm"
)
# Checking balance prior to matching
summary(m.out0)
sink("CEM_baseline.txt")
summary(m.out0, un = FALSE)
sink()  # returns output to the console

# Corsened exact matching
cutpoints <- list(
    Age = c(0, 20, 40, 60, 80, 100),
    TumorSize = "q3"
)
# Fit
m.out2 <- matchit(mymodel,
    data = df,
    method = "cem",
    estimand = "ATT", # Focus on ATE vs ATT
    cutpoints = cutpoints,
    ratio = 1
)
m.out2
sink("CEM_summary.txt")
summary(m.out2, un = FALSE)
sink()  # returns output to the console



png(
    file = "plots/check_matching1_CEM.png",
    units="in", width=5, height=5,
    res=300
)
plot(m.out2,
    type = "density", interactive = FALSE,
    which.xs = ~ Age + Gender + Race
)
dev.off()

png(
    file = "plots/check_matching2_CEM.png",
    units="in", width=5, height=5,
    res=300
)
plot(m.out2,
    type = "density", interactive = FALSE,
    which.xs = ~ TumorSize
)
dev.off()

png(
    file = "plots/check_matching3_CEM.png",
    units="in", width=5, height=5,
    res=300
)
plot(summary(m.out2))
dev.off()

m.data <- match.data(m.out2)
matched_data <- m.data

# Save data
write.csv(matched_data, "matched_data.csv", row.names = FALSE)

