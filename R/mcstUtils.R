#### Internal functions ####

  # .prepare_matrix
  # .get_irreducible
  # .get_cycle_complete
  # .get_cost
  # .get_arcs
  # .plot_mcstp


# Prepares the cost matrix from multiple input formats:
#   - Numeric vector
#   - Square numeric matrix or data.frame
#   - Edge list: data.frame or matrix with columns 'from', 'to',
#     and a cost column
#   - igraph object (weighted, undirected)
#   - Path to an Excel or CSV file (.xlsx/.xls or .csv): square
#     matrix (no headers) or edge list with columns 'from', 'to'

.prepare_matrix <- function(x) {

  # Numeric vector
  if (is.vector(x) && is.numeric(x)) {
    # Solve for N in the equation length(x) = N * (N - 1) / 2
    N <- (1 + sqrt(1 + 8 * length(x))) / 2

    if (abs(N - round(N)) > 1e-9)
      stop("Vector length does not correspond to N*(N-1)/2")

    N <- round(N)
    m <- matrix(0, N, N)
    m[lower.tri(m)] <- x
    m <- m + t(m)
    rownames(m) <- colnames(m) <- as.character(0:(N - 1L))

    return(m)
  }

  # Data frame or matrix
  if (is.data.frame(x) || is.matrix(x)) {
    x_df <- as.data.frame(x)
    cols <- tolower(colnames(x_df))

    is_edge <- ("from" %in% cols && "to" %in% cols)

    if (is_edge) { # Edge list format
      idx_f <- which(cols == "from")
      idx_t <- which(cols == "to")
      idx_c <- setdiff(seq_len(ncol(x_df)), c(idx_f, idx_t))[1]

      from  <- as.character(x_df[, idx_f])
      to    <- as.character(x_df[, idx_t])
      costs <- as.numeric(x_df[, idx_c])
      nodes <- sort(unique(c(from, to)))
      N     <- length(nodes)

      m <- matrix(0, N, N, dimnames = list(nodes, nodes))
      for (i in seq_along(from)) {
        m[from[i], to[i]] <- costs[i]
        m[to[i], from[i]] <- costs[i]
      }
      return(m)

    } else { # Square matrix format
      m <- as.matrix(x)
      storage.mode(m) <- "double"
      if (nrow(m) != ncol(m)) stop("Matrix must be square")
      rownames(m) <- colnames(m) <- as.character(0:(nrow(m) - 1L))

      return(m)
    }
  }

  # 'igraph' object
  if (inherits(x, "igraph")) {
    if (!requireNamespace("igraph", quietly = TRUE))
      stop("Package 'igraph' is required to handle igraph objects. Please install it")

    # Get adjacency matrix using weight attribute
    m <- as.matrix(igraph::as_adjacency_matrix(x, attr = "weight", sparse = FALSE))
    m <- pmax(m, t(m))
    rownames(m) <- colnames(m) <- as.character(seq_len(nrow(m)) - 1L)

    return(m)
  }

  # File path
  if (is.character(x) && length(x) == 1) {
    if (!file.exists(x))
      stop("File '", x, "' not found in the current working directory")

    ext <- tolower(tools::file_ext(x))

    if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE))
        stop("Package 'readxl' is required to read Excel files. Please install it")

      raw <- as.data.frame(readxl::read_excel(x, col_names = TRUE))

      # Check if first row is header or data
      if (!"from" %in% tolower(colnames(raw))) {
        raw <- as.data.frame(readxl::read_excel(x, col_names = FALSE))
        colnames(raw) <- NULL
      }

    } else if (ext == "csv") {

      line <- readLines(x, n = 1, warn = FALSE)
      df_func <- if (grepl(";", line)) read.csv2 else read.csv

      header_check <- df_func(x, header = FALSE, nrows = 1, stringsAsFactors = FALSE)
      has_head <- "from" %in% tolower(trimws(unlist(header_check)))

      raw <- df_func(x, header = has_head, stringsAsFactors = FALSE)
      if (!has_head) colnames(raw) <- NULL

    } else {
      stop("Unsupported file extension '", ext, "'. Accepted formats: .xlsx, .xls, .csv")
    }

    return(.prepare_matrix(raw))
  }

  stop("Unsupported input type")
}


# Computes the irreducible cost matrix barC from a cost matrix C
# The irreducible cost between node i and j is the maximum cost arc
# on the unique path connecting them in the MST
.get_irreducible <- function(C_mat) {

  N <- nrow(C_mat)
  C_irred <- C_mat

  # Finds the maximum arc on the optimal path between all pairs
  for (k in 1:N) {
    # Compares direct cost with the path via node k
    max_k <- outer(C_irred[, k], C_irred[k, ], pmax)
    C_irred <- pmin(C_irred, max_k)
  }
  diag(C_irred) <- 0

  return(C_irred)
}


# Computes the cycle-complete cost matrix C* from a cost matrix C
# The cycle-complete cost between i and j is the maximum irreducible
# cost across all subgraphs where a third agent k is removed
.get_cycle_complete <- function(C_mat) {

  N <- nrow(C_mat)
  C_cc <- matrix(0, N, N)
  dimnames(C_cc) <- dimnames(C_mat)

  # Precomputes irreducible matrices excluding each agent k
  sub_irreducibles <- list()
  agents <- 2:N

  for (k in agents) {
    sub_C <- C_mat[-k, -k]
    sub_irreducibles[[as.character(k)]] <- .get_irreducible(sub_C)
  }

  # Precompute full irreducible matrix for arc cases
  C_irred <- NULL

  # Computes only the upper triangle (C* is symmetric)
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {

      others <- setdiff(agents, c(i, j))

      if (length(others) == 0) {
        # Defaults to full irreducible cost if no other agents exist
        if (is.null(C_irred)) {
          C_irred <- .get_irreducible(C_mat)
        }
        max_val <- C_irred[i, j]
      } else {
        # Selects the maximum irreducible cost from the precomputed subgraphs
        max_val <- -Inf
        for (k in others) {
          sub_i <- if (i < k) i else i - 1
          sub_j <- if (j < k) j else j - 1
          max_val <- max(max_val, sub_irreducibles[[as.character(k)]][sub_i, sub_j])
        }
      }

      # Assign to both symmetrical positions
      C_cc[i, j] <- max_val
      C_cc[j, i] <- max_val
    }
  }

  return(C_cc)
}


# Optimized Prim's algorithm to compute the total cost of a MST
.get_cost <- function(sub_C) {

  N <- nrow(sub_C)
  # Return 0 if the coalition is empty or only contains the source
  if (N <= 1) return(0)

  # Initialize connection costs as infinite
  min_costs <- rep(Inf, N)
  min_costs[1] <- 0          # Starting point: the source
  visited <- rep(FALSE, N)   # Track connected nodes
  total <- 0

  for (p in 1:N) {

    # Select the cheapest unconnected node
    i <- which.min(ifelse(visited, Inf, min_costs))

    # Mark as connected and add its cost to the total
    visited[i] <- TRUE
    total <- total + min_costs[i]

    # Update best connection offers for all nodes using the new node 'i'
    # pmin compares the current best costs with the new available costs
    min_costs <- pmin(min_costs, sub_C[i, ])
  }

  return(total)
}


# Function that obtains the arcs of the MST for the full graph
# using Prim's algorithm to allow visualization
.get_arcs <- function(C_mat) {
  # n <- nrow(C_mat) - 1
  S <- "0"
  # N_minus_S <- as.character(1:n)
  N_minus_S <- setdiff(rownames(C_mat), "0")
  n <- length(N_minus_S)

  arcs <- data.frame(
    i = character(n),
    j = character(n),
    stage = integer(n),
    stringsAsFactors = FALSE
  )

  for (p in 1:n) {
    sub_C <- C_mat[S, N_minus_S, drop = FALSE]
    min_cost <- min(sub_C)
    res <- which(sub_C == min_cost, arr.ind = TRUE)

    i <- S[res[1, 1]]
    j <- N_minus_S[res[1, 2]]

    arcs[p, ] <- list(i, j, p)

    S <- c(S, j)
    N_minus_S <- setdiff(N_minus_S, j)
  }
  return(arcs)
}


# General plotting function for mcstp problems
.plot_mcstp <- function(adj, mst_arcs, allocation, main_title, sub_title, show_stage = TRUE) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' required for drawing. Please install it")
  }

  n <- nrow(adj) - 1
  adj_plot <- adj
  adj_plot[is.infinite(adj_plot)] <- 0
  g <- igraph::graph_from_adjacency_matrix(adj_plot,
                                           mode = "undirected", weighted = TRUE)

  arc_colors <- rep("grey80", igraph::ecount(g))
  arc_widths <- rep(1, igraph::ecount(g))
  arc_labels <- as.character(igraph::E(g)$weight)

  arcs <- igraph::as_edgelist(g, names = TRUE)

  for(k in 1:nrow(mst_arcs)) {
    i <- as.character(mst_arcs$i[k])
    j <- as.character(mst_arcs$j[k])
    stage <- NULL
    if (show_stage && "stage" %in% colnames(mst_arcs)) {
      stage <- mst_arcs$stage[k]
    }

    for(l in 1:nrow(arcs)) {
      if((arcs[l,1] == i && arcs[l,2] == j) ||
         (arcs[l,1] == j && arcs[l,2] == i)) {
        arc_colors[l] <- "red"
        arc_widths[l] <- 2
        if (!is.null(stage)) {
          arc_labels[l] <- paste0(arc_labels[l], " [", stage, "]")
        }
      }
    }
  }

  v_labels <- as.character(0:n)
  v_colors <- c("#F0FFFF", rep("#D8E6E6", n))

  for(i in 1:n) {
    if (is.null(allocation)) {
      v_labels[i+1] <- as.character(i)
    } else {
      v_labels[i+1] <- paste0(i, "\n(", round(allocation[i], 2), ")")
    }
  }

  plot(g,
       layout = igraph::layout_in_circle(g),
       edge.label = arc_labels,
       edge.label.color = "#27408B",
       edge.label.cex = 0.9,
       edge.color = arc_colors,
       edge.width = arc_widths,
       edge.lty = 1,
       vertex.label = v_labels,
       vertex.color = v_colors,
       vertex.frame.color = "black",
       vertex.label.color = "black",
       vertex.label.font = 1,
       vertex.label.cex = 0.8,
       vertex.size = 40,
       main = if(main_title == "") NULL else main_title)

  if(sub_title != "") {
    mtext(sub_title, side = 3, line = 0.5, cex = 0.7,
          family = "sans", font = 3, col = "#444444")
  }
}

