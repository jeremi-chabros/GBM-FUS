library(MatchIt)
library("marginaleffects")
library(survival)
library(lubridate)

# Read data
read_csv <- read.csv("new.csv")
df <- na.omit(read_csv) # Drop missing values
df <- df[!(is.na(df$DeathCensorDate) | df$DeathCensorDate == ""), ]
df <- df[!(df$X1p19q == 1), ]
df <- df[df$Chemotherapy == 1, ]
df <- df[df$Radiotherapy == 1, ]
df <- df[!(df$TumorSize == "Missing"), ]

# Dichotomize race
df$Race <- ifelse(df$Race == "White", "White", "Non-white")
# mymodel <- as.formula(FUS ~ Age + Gender + Race + MGMT + IDH + X1p19q)
mymodel <- as.formula(FUS ~ Age + Gender + Race + TumorSize)

m.out0 <- matchit(mymodel,
    data = df,
    method = NULL, distance = "glm"
)
sink("optimal_baseline.txt")
summary(m.out0)
sink()

# # Full matching on a probit PS
# m.out_optimal <- matchit(mymodel,
#     data = df,
#     method = "full",
#     distance = "glm",
#     ratio = 1,
#     replace = FALSE,
#     link = "probit"
# )

# Optimal matching
m.out_optimal <- matchit(mymodel,
                         data = df,
                         method = "optimal",
                         distance = "glm",
                         ratio = 1)

# Check balance to assess the quality of the matches
sink("optimal_summary.txt")
summary(m.out_optimal)
sink()


png(
    file = "plots/check_matching1_optimal.png",
    units="in", width=5, height=5,
    res=300
)
plot(m.out_optimal,
    type = "density", interactive = FALSE,
    which.xs = ~ Age + Gender + Race
)
dev.off()

png(
    file = "plots/check_matching2_optimal.png",
    units="in", width=5, height=5,
    res=300
)
plot(m.out_optimal,
    type = "density", interactive = FALSE,
    which.xs = ~ TumorSize
)
dev.off()

png(
    file = "plots/check_matching3_optimal.png",
    units="in", width=5, height=5,
    res=300
)
plot(summary(m.out_optimal))
dev.off()

# Extract the matched data and save only the matched subset to a new CSV file
m.data_optimal <- match.data(m.out_optimal)
write.csv(m.data_optimal, "matched_optimal.csv", row.names = FALSE)


