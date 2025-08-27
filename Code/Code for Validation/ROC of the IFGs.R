# Load necessary libraries
library(pROC)
library(ggplot2)

# Load the dataset
file_path <- "F:\\Green University\\Personal\\MSc in RUET\\MSc Thesis\\Thesis Work\\Journal\\Validation\\Zscores_and_stages_of_influential_genes.csv"
dataset <- read.csv(file_path)

# Remove Patient_ID column
dataset <- dataset[, -1]

# Convert Stage to binary (0 = Early, 1 = Late)
dataset$Stage <- as.numeric(as.factor(dataset$Stage)) - 1  

# Initialize an empty list to store ROC curves
roc_curves <- list()

# Loop through all 9 genes and compute ROC curves
genes <- colnames(dataset)[-1]  # Exclude "Stage"

for (gene in genes) {
  roc_curve <- roc(dataset$Stage, dataset[[gene]])
  roc_smooth <- smooth(roc_curve)  # Smooth ROC curve
  roc_curves[[gene]] <- roc_smooth
}

# Plot all ROC curves
plot(roc_curves[[1]], col = 1, lwd = 2, main = "Smoothed ROC Curves for All Genes")

# Add remaining ROC curves
for (i in 2:length(genes)) {
  lines(roc_curves[[i]], col = i, lwd = 2)
}

# Add legend
legend("bottomright", legend = genes, col = 1:length(genes), lwd = 2)
