setwd("C:/Users/Quinc/OneDrive - University of Kentucky/R/mcc_method_comparison/code")
library(mapproj)
library(terra)
library(httr)
library(tidyr)
library(dplyr)
library(reshape2)


testing_scenario <- "clustering_mixtureautofeatures_seeds_1"
path <- paste("../results/", testing_scenario, sep = "")

# List all .rds files in the folder
files <- list.files(path, pattern = "\\.rds$")

# Read all .rds files into a list
data_list <- lapply(paste0(path, "/", files), readRDS)

# Initialize an empty dataframe
combined_df <- data.frame()

# Loop over each file
for(i in c(1:40, 121:140, 221:240, 41:80, 181:220, 81:120, 141:180)) {
  data.name <- data_list[[i]]$data.name
  clustering.method <-data_list[[i]]$clustering.method
  feature.selection.method <- as.character(data_list[[i]]$feature.selection.method)
  ARI <- data_list[[i]]$ARI
  Purity <- data_list[[i]]$Purity
  Jaccard <- data_list[[i]]$Jaccard
  FM <- data_list[[i]]$FM
  
  # Combine these into a dataframe
  df <- data.frame( data.name, clustering.method, feature.selection.method, ARI, Purity, Jaccard, FM)
  
  # Add this dataframe to the combined dataframe
  combined_df <- rbind(combined_df, df)
  rownames(combined_df) = NULL
}

# Now, `combined_df` is a dataframe that contains all the data from your RDS files

# Load the necessary library
library(ggplot2)

# Assuming `combined_df` is your dataframe and it has the same structure as the data in the image
# Convert factors to character for proper ordering in the plot
combined_df$data.name <- factor(as.character(combined_df$data.name), levels = unique(combined_df$data.name))
combined_df$feature.selection.method  <- gsub("gm", "DIFS", combined_df$feature.selection.method)
  
combined_df$feature.selection.method <- as.character(combined_df$feature.selection.method)

combined_df$feature.selection.method <- factor(combined_df$feature.selection.method, levels = c("SC3", "Seurat", "Feast", "SCT", "DIFS"))

# Create the plot
# Determine the number of unique feature selection methods
num_methods <- length(unique(combined_df$feature.selection.method))


library(ggplot2)



# Plot
p <- ggplot(combined_df, aes(x=clustering.method, y=ARI, fill=feature.selection.method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~data.name, ncol=4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x="Clustering methods", y="Adjusted Rand Index (ARI)") +
  scale_fill_brewer(palette="Spectral")

print(p)


ggsave(paste("../results/plotout/", testing_scenario, ".pdf", sep = ""), plot = p, height = 8, width = 8)



combined_df_ari <- combined_df[, -(5:7)]
combined_df_ari$method_combined <- with(combined_df_ari, paste(clustering.method, feature.selection.method, sep = "_"))

# Calculate average ARI for each method combination
average_ari_methods <- combined_df_ari %>%
  group_by(method_combined) %>%
  summarise(Avg_ARI_Method = mean(ARI, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(Avg_ARI_Method))

# Calculate average ARI for each dataset
average_ari_datasets <- combined_df_ari %>%
  group_by(data.name) %>%
  summarise(Avg_ARI_Dataset = mean(ARI, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Avg_ARI_Dataset)

# Update factor levels based on average ARI
combined_df_ari$method_combined <- factor(combined_df_ari$method_combined, levels = average_ari_methods$method_combined)
combined_df_ari$data.name <- factor(combined_df_ari$data.name, levels = average_ari_datasets$data.name)

# Melt 'combined_df' into a long format suitable for plotting
data_long <- melt(combined_df_ari[, -(2:3)], id.vars = c("data.name", "method_combined"), variable.name = "Metric", value.name = "ARI")


library(viridis)  # Load the viridis package for vibrant color palettes



# Correcting the function call for viridis in the heatmap
heatmap_plot <- ggplot(data_long, aes(x = method_combined, y = data.name, fill = ARI)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "viridis") +  # Corrected usage of viridis palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(x = "Method Combination", fill = "ARI")

# Display the corrected heatmap
print(heatmap_plot)

# Save the heatmap with viridis colors
ggsave("../results/plotout/heatmap_viridis.pdf", plot = heatmap_plot, height = 6, width = 10)





library(ggplot2)
library(viridis)  # Load the viridis package for vibrant color palettes
library(dplyr)
library(reshape2)
library(patchwork)  # Library for combining plots

unique_methods <- unique(combined_df_ari$clustering.method)

# Initialize a list to store plots
plot_list <- list()
plot_count <- 1

for(method in unique_methods[1:2]) {
  # Filter data for the current method
  filtered_data <- combined_df_ari %>%
    filter(clustering.method == method)
  
  # Calculate average ARI for each method combination
  average_ari_methods_filtered <- filtered_data %>%
    group_by(feature.selection.method) %>%
    summarise(Avg_ARI_Method = mean(ARI, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(Avg_ARI_Method))
  
  # Calculate average ARI for each dataset
  average_ari_datasets_filtered <- combined_df_ari %>%
    group_by(data.name) %>%
    summarise(Avg_ARI_Dataset = mean(ARI, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(Avg_ARI_Dataset)
  
  # Ensure the factor levels are set correctly within the loop for filtered_data
  filtered_data$feature.selection.method <- factor(filtered_data$feature.selection.method, levels = average_ari_methods_filtered$feature.selection.method)
  filtered_data$data.name <- factor(filtered_data$data.name, levels = average_ari_datasets_filtered$data.name)
  
  # Melt filtered_data into a long format suitable for plotting
  data_long <- melt(filtered_data[, -c(2, 5)], id.vars = c("data.name", "feature.selection.method"), variable.name = "Metric", value.name = "ARI")
  
  # Create a heatmap plot using the viridis color palette with adjusted legend size
  p <- ggplot(data_long, aes(x = data.name, y = feature.selection.method, fill = ARI)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c() +  # Use viridis palette
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1)),  # Reduce font size of x-axis labels
      axis.text.y = element_text(size = rel(1)),  # Reduce font size of y-axis labels
      axis.title.x = element_text(size = rel(1)),  # Reduce font size of x-axis title
      axis.title.y = element_text(size = rel(1)),  # Reduce font size of y-axis title
      legend.text = element_text(size = rel(0.3)),  # Further reduce font size in the legend
      legend.title = element_text(size = rel(0.3)),  # Reduce font size of the legend title
      plot.title = element_text(size = rel(1))  # Reduce font size of the plot title
    ) +
    labs(title = method, x = "Dataset", y = "Feature Selection Method")
  
  # Add the plot to the list
  plot_list[[plot_count]] <- p + ylab(NULL) + theme(
     legend.position = "none",
     axis.text.x = element_blank(),    # Remove x-axis tick labels
     axis.ticks.x = element_blank()    # Remove x-axis ticks
  ) +xlab(NULL)
  plot_count <- plot_count + 1
}

for(method in unique_methods[3]) {
  # Filter data for the current method
  filtered_data <- combined_df_ari %>%
    filter(clustering.method == method)
  
  # Calculate average ARI for each method combination
  average_ari_methods_filtered <- filtered_data %>%
    group_by(feature.selection.method) %>%
    summarise(Avg_ARI_Method = mean(ARI, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(Avg_ARI_Method))
  
  # Calculate average ARI for each dataset
  average_ari_datasets_filtered <- combined_df_ari %>%
    group_by(data.name) %>%
    summarise(Avg_ARI_Dataset = mean(ARI, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(Avg_ARI_Dataset)
  
  # Ensure the factor levels are set correctly within the loop for filtered_data
  filtered_data$feature.selection.method <- factor(filtered_data$feature.selection.method, levels = average_ari_methods_filtered$feature.selection.method)
  filtered_data$data.name <- factor(filtered_data$data.name, levels = average_ari_datasets_filtered$data.name)
  
  # Melt filtered_data into a long format suitable for plotting
  data_long <- melt(filtered_data[, -c(2, 5)], id.vars = c("data.name", "feature.selection.method"), variable.name = "Metric", value.name = "ARI")
  
  # Create a heatmap plot using the viridis color palette with adjusted legend size
  p <- ggplot(data_long, aes(x = data.name, y = feature.selection.method, fill = ARI)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c() +  # Use viridis palette
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1)),  # Reduce font size of x-axis labels
      axis.text.y = element_text(size = rel(1)),  # Reduce font size of y-axis labels
      axis.title.x = element_text(size = rel(1)),  # Reduce font size of x-axis title
      axis.title.y = element_text(size = rel(1)),  # Reduce font size of y-axis title
      legend.text = element_text(size = rel(0.3)),  # Further reduce font size in the legend
      legend.title = element_text(size = rel(0.3)),  # Reduce font size of the legend title
      plot.title = element_text(size = rel(1))  # Reduce font size of the plot title
    ) +
    labs(title = method, x = "Dataset", y = "Feature Selection Method")
  
  # Add the plot to the list
  plot_list[[plot_count]] <- p + ylab(NULL) + theme(
    #legend.position = "none",
    axis.text.x = element_blank(),    # Remove x-axis tick labels
    axis.ticks.x = element_blank()    # Remove x-axis ticks
  ) +xlab(NULL)
  plot_count <- plot_count + 1
}


for(method in unique_methods[4]) {
  # Filter data for the current method
  filtered_data <- combined_df_ari %>%
    filter(clustering.method == method)
  
  # Calculate average ARI for each method combination
  average_ari_methods_filtered <- filtered_data %>%
    group_by(feature.selection.method) %>%
    summarise(Avg_ARI_Method = mean(ARI, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(Avg_ARI_Method))
  
  # Calculate average ARI for each dataset
  average_ari_datasets_filtered <- combined_df_ari %>%
    group_by(data.name) %>%
    summarise(Avg_ARI_Dataset = mean(ARI, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(Avg_ARI_Dataset)
  
  # Ensure the factor levels are set correctly within the loop for filtered_data
  filtered_data$feature.selection.method <- factor(filtered_data$feature.selection.method, levels = average_ari_methods_filtered$feature.selection.method)
  filtered_data$data.name <- factor(filtered_data$data.name, levels = average_ari_datasets_filtered$data.name)
  
  data_long <- melt(filtered_data[, -c(2, 5)], id.vars = c("data.name", "feature.selection.method"), variable.name = "Metric", value.name = "ARI")
  
  # Create a heatmap plot using the viridis color palette with adjusted legend size
  p <- ggplot(data_long, aes(x = data.name, y = feature.selection.method, fill = ARI)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c() +  # Use viridis palette
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1)),  # Reduce font size of x-axis labels
      axis.text.y = element_text(size = rel(1)),  # Reduce font size of y-axis labels
      axis.title.x = element_text(size = rel(1)),  # Reduce font size of x-axis title
      axis.title.y = element_text(size = rel(1)),  # Reduce font size of y-axis title
      legend.text = element_text(size = rel(0.3)),  # Further reduce font size in the legend
      legend.title = element_text(size = rel(0.3)),  # Reduce font size of the legend title
      plot.title = element_text(size = rel(1))  # Reduce font size of the plot title
    ) +
    labs(title = method, x = "Dataset", y = "Feature Selection Method")
  
  # Add the plot to the list
  plot_list[[plot_count]] <- p + ylab(NULL) + theme(
     legend.position = "none",
     axis.text.x = element_blank(),    # Remove x-axis tick labels
     axis.ticks.x = element_blank()    # Remove x-axis ticks
  ) +xlab(NULL)
  plot_count <- plot_count + 1
}

for(method in unique_methods[5]) {
  # Filter data for the current method
  filtered_data <- combined_df_ari %>%
    filter(clustering.method == method)
  
  # Calculate average ARI for each method combination
  average_ari_methods_filtered <- filtered_data %>%
    group_by(method_combined) %>%
    summarise(Avg_ARI_Method = mean(ARI, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(Avg_ARI_Method))
  
  # Calculate average ARI for each dataset
  average_ari_datasets_filtered <- combined_df_ari %>%
    group_by(data.name) %>%
    summarise(Avg_ARI_Dataset = mean(ARI, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(Avg_ARI_Dataset)
  
  # Calculate average ARI for each method combination
  average_ari_methods_filtered <- filtered_data %>%
    group_by(feature.selection.method) %>%
    summarise(Avg_ARI_Method = mean(ARI, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(Avg_ARI_Method))
  
  # Calculate average ARI for each dataset
  average_ari_datasets_filtered <- combined_df_ari %>%
    group_by(data.name) %>%
    summarise(Avg_ARI_Dataset = mean(ARI, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(Avg_ARI_Dataset)
  
  # Ensure the factor levels are set correctly within the loop for filtered_data
  filtered_data$feature.selection.method <- factor(filtered_data$feature.selection.method, levels = average_ari_methods_filtered$feature.selection.method)
  filtered_data$data.name <- factor(filtered_data$data.name, levels = average_ari_datasets_filtered$data.name)
  
  data_long <- melt(filtered_data[, -c(2, 5)], id.vars = c("data.name", "feature.selection.method"), variable.name = "Metric", value.name = "ARI")
  
  # Create a heatmap plot using the viridis color palette with adjusted legend size
  p <- ggplot(data_long, aes(x = data.name, y = feature.selection.method, fill = ARI)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c() +  # Use viridis palette
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1)),  # Reduce font size of x-axis labels
      axis.text.y = element_text(size = rel(1)),  # Reduce font size of y-axis labels
      axis.title.x = element_text(size = rel(1)),  # Reduce font size of x-axis title
      axis.title.y = element_text(size = rel(1)),  # Reduce font size of y-axis title
      legend.text = element_text(size = rel(0.3)),  # Further reduce font size in the legend
      legend.title = element_text(size = rel(0.3)),  # Reduce font size of the legend title
      plot.title = element_text(size = rel(1))  # Reduce font size of the plot title
    ) +
    labs(title = method, x = "Dataset", y = "Feature Selection Method")
  
  # Add the plot to the list
  plot_list[[plot_count]] <- p + ylab(NULL) + theme(
    legend.position = "none",
    #axis.text.x = element_blank(),    # Remove x-axis tick labels
    #axis.ticks.x = element_blank()    # Remove x-axis ticks
  ) #+xlab(NULL)
  plot_count <- plot_count + 1
}

# Combine all plots into one using the patchwork package
combined_plot <- wrap_plots(plot_list, ncol = 1) + plot_layout(guides = "collect")  # Collect all legends into one

# Display the combined plot with a single, collected legend
print(combined_plot)

# Save the combined heatmap plot with a single, smaller legend
ggsave("../results/plotout/combined_named_heatmap_small_font.pdf", plot = combined_plot, height = 2 * length(unique_methods), width = 5)






library(dplyr)
library(ggplot2)
library(tidyr)

# Ensure combined_df is defined or loaded before proceeding
# Define high contrast colors for the metrics
high_contrast_colors <- c("Avg_ARI" = "#007FFF", "Avg_FM" = "#FF0038", "Avg_Purity" = "#50C878", "Avg_Jaccard" = "#FFBF00")

# Combine clustering and feature selection methods into a new column
combined_df$method_combined <- with(combined_df, paste(clustering.method, feature.selection.method, sep = "_"))

# Compute averages for metrics grouped by feature selection method
averages_df <- combined_df %>%
  group_by(feature.selection.method) %>%
  summarise(
    Avg_ARI = mean(ARI, na.rm = TRUE),
    Avg_FM = mean(FM, na.rm = TRUE),
    Avg_Purity = mean(Purity, na.rm = TRUE),
    Avg_Jaccard = mean(Jaccard, na.rm = TRUE),
    .groups = 'drop'  # This ensures the summarise function does not complain about row information
  )

# Pivot data for plotting
long_df <- pivot_longer(averages_df, 
                        cols = starts_with("Avg_"), 
                        names_to = "Metric", 
                        values_to = "Value")

# Create a line chart with the specified high-contrast colors
line_chart <- ggplot(long_df, aes(x = feature.selection.method, y = Value, color = Metric, group = Metric)) +
  geom_line(size = 1) +  # Define line size
  geom_point(size = 4) +  # Define point size
  scale_color_manual(values = high_contrast_colors) +  # Apply custom colors
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Feature Selection Method", y = "Average Metric Value", title = "Model Performance Across Different Metrics")

line_chart  # Display the plot


# Save the plot with a smaller size
output_dir <- "../results/plotout/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
ggsave(paste0(output_dir, "line_chart.pdf"), plot = line_chart, height = 5, width = 7)  # Adjusted size

# Function to create and save a line plot
create_line_plot <- function(df, scenario_name) {
  # Compute averages for metrics grouped by clustering method
  averages_df <- df %>%
    group_by(clustering_method) %>%
    summarise(
      Avg_ARI = mean(ARI, na.rm = TRUE),
      Avg_FM = mean(FM, na.rm = TRUE),
      Avg_Purity = mean(Purity, na.rm = TRUE),
      Avg_Jaccard = mean(Jaccard, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Pivot data for plotting
  long_df <- pivot_longer(averages_df, 
                          cols = starts_with("Avg_"), 
                          names_to = "Metric", 
                          values_to = "Value")
  
  # Create a line chart with the specified high-contrast colors
  line_chart <- ggplot(long_df, aes(x = clustering_method, y = Value, color = Metric, group = Metric)) +
    geom_line(size = 1) +  # Define line size
    geom_point(size = 4) +  # Define point size
    scale_color_manual(values = high_contrast_colors) +  # Apply custom colors
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Clustering Method", y = "Average Metric Value", title = paste("Model Performance Across Different Metrics -", scenario_name))
  
  # Save the plot
  output_dir <- "../results/plot/"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(paste0(output_dir, "line_chart_", scenario_name, ".pdf"), plot = line_chart, height = 5, width = 7)
  
  # Print the line chart
  print(line_chart)
}

# Create and save line plots for each scenario
create_line_plot(realdata_df, "realdata")
create_line_plot(heterogeneous_df, "heterogeneous")
create_line_plot(homogeneous_df, "homogeneous")

# Print a message indicating completion
print("Line plots saved for each scenario")


