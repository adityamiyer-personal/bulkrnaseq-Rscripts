#https://x.com/i/grok/share/Nr32mi5X26id5NrcDfGtI4VaL
#check by running the script
# Assuming 'expr' is your normalized RNA-seq count matrix and 'group' is your group vector

library(ggplot2)

# Step 1: Compute MAD scores
group_levels <- unique(group)
median_profiles <- lapply(group_levels, function(g) {
  samples_g <- which(group == g)
  apply(expr[, samples_g, drop = FALSE], 1, median)
})
names(median_profiles) <- group_levels

mad_scores <- numeric(ncol(expr))
for (j in 1:ncol(expr)) {
  g <- group[j]
  median_profile <- median_profiles[[g]]
  abs_dev <- abs(expr[, j] - median_profile)
  mad_scores[j] <- median(abs_dev)
}

# Step 2: Identify outliers
outliers <- logical(ncol(expr))
for (g in group_levels) {
  samples_g <- which(group == g)
  mad_g <- mad_scores[samples_g]
  median_mad <- median(mad_g)
  mad_mad <- median(abs(mad_g - median_mad))
  threshold <- median_mad + 3 * mad_mad
  outliers[samples_g] <- mad_g > threshold
}

# Create data frame
plot_data <- data.frame(
  sample = colnames(expr),
  group = group,
  mad_score = mad_scores,
  outlier = outliers
)

# Create scatter plot
ggplot(plot_data, aes(x = group, y = mad_score, color = outlier)) +
  geom_point(position = position_jitter(width = 0.2)) +
  theme_minimal() +
  labs(title = "MAD Scores by Sample Group", x = "Group", y = "MAD Score") +
  scale_color_manual(values = c("black", "red"), labels = c("Normal", "Outlier")) +
  guides(color = guide_legend(title = "Status"))
