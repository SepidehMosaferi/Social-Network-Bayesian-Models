# Network diagrams
library(igraph)

network_folder <- "/Users/home/Library/CloudStorage/Application"
selected_village <- 15
input_file <- file.path(network_folder,paste0("Vil15.RData")) #change this
output_file <- file.path(network_folder,
                          paste0("village", selected_village, "-networks-no-isolates.jpeg"))

if (!file.exists(input_file)) stop("File not found: ", input_file)
env <- new.env(parent = emptyenv())
load(input_file, envir = env)

relation_names <- c("nonrel", "rel", "templecompany")
matrices <- lapply(relation_names, function(nm) {
  key <- paste0("Matrix_", nm)
  if (!exists(key, envir = env, inherits = FALSE)) stop("Missing matrix: ", key)
  get(key, envir = env, inherits = FALSE)
})
names(matrices) <- relation_names

make_graph <- function(mat, relation) {
  if (!is.matrix(mat) || nrow(mat) != ncol(mat)) {
    stop(relation, ": matrix must be square")
  }
  if (is.null(rownames(mat)) || is.null(colnames(mat)) ||
      !identical(rownames(mat), colnames(mat))) {
    stop(relation, ": matrix row and column IDs must match in order")
  }
  # For visualization only
  adjacency <- (mat == 1)
  adjacency[is.na(adjacency)] <- FALSE
  diag(adjacency) <- FALSE
  g <- graph_from_adjacency_matrix(adjacency, mode = "directed", diag = FALSE)

  # Remove only nodes with BOTH indegree=0 AND outdegree=0.
  isolated <- degree(g, mode = "in") == 0 & degree(g, mode = "out") == 0
  cat(sprintf("%s: %d actors; %d complete isolates removed; %d plotted; %d directed edges\n",
              relation, vcount(g), sum(isolated), sum(!isolated), ecount(g)))
  delete_vertices(g, V(g)[isolated])
}

graphs <- Map(make_graph, matrices, relation_names)

# Save figure
png(output_file, width = 1800, height = 750, res = 150, bg = "white")
oldpar <- par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))
set.seed(1234)
for (relation in relation_names) {
  g <- graphs[[relation]]
  if (vcount(g) == 0) {
    plot.new()
    title(main = sprintf("%s (Village %s)\n(no observed ties)",
                         relation, selected_village))
    next
  }
  # Fruchterman-Reingold layout (independent layout per relationship).
  xy <- layout_with_fr(as.undirected(g, mode = "collapse"))
  plot(g,
       layout = xy,
       main = sprintf("%s (Village %s)", relation, selected_village),
       vertex.size = 3.8,
       vertex.color = "#e75a76",
       vertex.frame.color = "black",
       vertex.label = NA,
       edge.color = "black",
       edge.width = 0.65,
       edge.arrow.size = 0.18,
       edge.curved = 0,
       margin = 0.05)
}
par(oldpar)
dev.off()
cat("Saved figure: ", output_file, "\n", sep = "")
