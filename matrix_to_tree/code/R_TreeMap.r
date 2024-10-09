 # install.packages("ape")

list.files()

# Define a function that accepts the number of clusters k and plots the dendrogram
plot_dendrogram_with_colors <- function(
    dm, 
    k=4,
    method = "complete",
    type="phylogram",
    cex = 1, 
    label_offset = 0.02
) {
  # Convert dm to a dist object if it is not already
  if (!inherits(dm, "dist")) {
    dm_dist <- as.dist(as.matrix(dm))
  } else {
    dm_dist <- dm
  }

  # Compute the hierarchical clustering using the dist object
  hc <- hclust(dm_dist, method = method)

  # Determine cluster membership for each item
  clusters <- cutree(hc, k)

  # Generate k distinct colors
  colors <- rainbow(k)

  # Map the clusters to colors
  cluster_colors <- colors[clusters]

  # Load the 'ape' package; install it if necessary with install.packages("ape")
  library(ape)
    
  # Convert hc to class phylo
  phylo_tree <- as.phylo(hc)
 
  # Plot the tree with type "radial"
  if (type == "radial") {
    par(mar=c(1,1,1,1)+0.1)  # Adjust the margin to make sure there is space for labels
    plot(phylo_tree, type=type, tip.color=cluster_colors, cex = cex, label.offset = label_offset)

    # Optionally, you can add another plotting line to display the labels
    # with more control and customization, e.g.:
    add.phylo.labels(phylo_tree, all=TRUE, cex=cex, align=TRUE)
  } else {
    # For other plot types
    plot(phylo_tree, type=type, tip.color=cluster_colors, cex = cex, label.offset = label_offset)
  }       
}

list.files()

# Read the input data from a CSV file
similarity_matrix <- read.csv(file="FtsZ_inhibitor.csv", header=TRUE, sep=",", row.names=1)

# Scaling the similarity matrix to the range of 0 to 1
similarity_matrix_scaled <- similarity_matrix / 100

# Computing the distance matrix by subtracting the scaled similarity from 1
dm <- 1 - similarity_matrix_scaled

# Now call the function with the desired number of clusters
desired_k <- 4 # Replace 4 with your desired number of clusters

plot_dendrogram_with_colors(
    dm, 
    desired_k,
    method = "complete",
    type="phylogram"    
)

plot_dendrogram_with_colors(
    dm, desired_k,
    method = "complete",
    type="cladogram"    
)

plot_dendrogram_with_colors(
    dm, desired_k,
    method = "complete",
    type="fan"    
)

# plot_dendrogram_with_colors(
#     dm, desired_k,
#     method = "complete",
#     type="radial"    
# )

# Read the input data from a CSV file
similarity_matrix <- read.csv(file="matrix_of_jinglan_v2.csv", header=TRUE, sep=",", row.names=1)

# Scaling the similarity matrix to the range of 0 to 1
similarity_matrix_scaled <- similarity_matrix / 100

# Computing the distance matrix by subtracting the scaled similarity from 1
dm <- 1 - similarity_matrix_scaled

# Now call the function with the desired number of clusters
desired_k <- 4 # Replace 4 with your desired number of clusters

plot_dendrogram_with_colors(
    dm, 
    desired_k,
    method = "complete",
    type="phylogram"    
)
plot_dendrogram_with_colors(
    dm, desired_k,
    method = "complete",
    type="cladogram"    
)

plot_dendrogram_with_colors(
    dm, desired_k,
    method = "complete",
    type="fan"    
)


# plot_dendrogram_with_colors(
#     dm, desired_k,
#     method = "complete",
#     type="radial"    
# )
