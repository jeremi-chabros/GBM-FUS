# =========================================================================
# Author: Jeremi Chabros
# Affiliations:
# 1) University of Cambridge School of Clinical Medicine
# 2) Computational Neuroscience Outcomes Center, Department of Neurosurgery,
#    Brigham & Women’s Hospital, Harvard Medical School
# Date: 09/09/2023
# Description: This script prepares data by matching cases using
# Coarsened Exact Matching in the MatchIt package, then saves the matched
# data, summary statistics, and density plots.
# To avoid errors, it is recommended to run line-by-line or run twice as whole.
# =========================================================================

# Load required libraries
{
    library(MatchIt)
    library("marginaleffects")
    library(survival)
    library(lubridate)
    library(ggplot2)
    library("ggpubr")
    library(gridExtra)
    library(mice)
    library(Hmisc)
    library(dplyr)
    library(cobalt)
}
# Set relative directory for the project
# setwd("path/to/your/directory")
setwd("/Users/jjc/Research/GBM-FUS/")

# =========================================================================
# Functions
# =========================================================================

# Define functions to read and preprocess data
read_and_process_data <- function(file_path) {
    df <- read.csv(file_path)
    df <- subset(df, !(is.na(DeathCensorDate) | DeathCensorDate == "Unknown"))
    df <- subset(df, Chemotherapy == 1)
    df <- subset(df, Radiotherapy == 1)
    df <- subset(df, KPSgeq70 == 1)
    df <- subset(df, Age <= 80)
    df$MGMT[df$MGMT == "Unknown" & df$FUS == 1] <- "plchldr"
    # df <- subset(df, MGMT != "plchldr")
    df <- subset(df, Location %in% c("Frontal", "Temporal", "Parietal", "Occipital"))
    df <- subset(df, MGMT != "Unknown")
    df$Race <- as.factor(ifelse(df$Race == "White", "White", "Non-white")) # Dichotomize race
    df$Age <- as.numeric(df$Age)
    df$Gender <- as.factor(df$Gender)
    df$MGMT <- as.factor(df$MGMT)
    df$IDH <- as.factor(df$IDH)
    df$Location <- as.factor(df$Location)
    df$Dead <- as.factor(df$Dead)
    df$Progression <- as.factor(df$Progression)
    df$FUS <- as.factor(df$FUS)
    df$TumorSize <- as.numeric(df$TumorSize)
    mode_MGMT <- names(sort(table(df$MGMT), decreasing = TRUE))[1]
    df$MGMT[df$MGMT == "plchldr"] <- "Methylated"
    return(df)
}

# Define function for matching
perform_matching <- function(df, formula, method, distance, cutpoints) {
    m.out <- matchit(formula,
        data = df,
        method = method,
        distance = distance,
        cutpoints = cutpoints
    )
    return(m.out)
}

# Function to save density plots
save_density_plot <- function(m.out, which_vars, file_path) {
    png(file = file_path, units = "in", width = 5, height = 5, res = 300)
    plot(m.out, type = "density", interactive = FALSE, which.xs = which_vars)
    title(main = "")
    dev.off()
}

# =========================================================================
# Main code starts here
# =========================================================================

# Read and preprocess data
file_path <- "data/FinalData.csv"
df <- read.csv(file_path)
df <- read_and_process_data(file_path)

# Define model formula
mymodel <- as.formula(FUS ~ Age + IDH + MGMT)

# Check balance prior to matching
m.out0 <- matchit(as.formula(FUS ~ Age + Gender + Race + IDH + MGMT + TumorSize),
    data = df,
    method = NULL, distance = "glm"
)
m.out0
summary(m.out0)

# Perform matching
cutpoints <- list(
    Age = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
)
m.out2 <- perform_matching(df, mymodel, "cem", "glm", cutpoints)
m.out2
summary(m.out2)
matched_data <- match.data(m.out2)
bal.tab(m.out2, data = m.out2, all.vars = TRUE, un = TRUE) # print balance measures table (unadjusted vs. adjusted)

# Save summary to a text file
sink("results/CEM_summary.txt")
summary(m.out2, un = FALSE)
sink()

# Plot and save matching results
save_density_plot(m.out2, ~ Age + MGMT + IDH, "results/matching_check.png")

# Save matched data to a CSV file
matched_data <- match.data(m.out2)
write.csv(matched_data, "data/MatchedData.csv", row.names = FALSE)

# =========================================================================