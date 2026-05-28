# RULES DEFINED THROUGH THE PROBLEM

  # Based on Prim's algorithm
    # Bird rule:                          mcstBird
    # Dutta-Kar rule:                     mcstDuttaKar

  # Based on Kruskal's algorithm
    # Folk rule:                          mcstFolk
    # Optimistic weighted Shapley rule:   mcstOWShapley
    # Pessimistic weighted Shapley rule:  mcstPWShapley

  # Based on Boruvka's algorithm
    # Boruvka rule:                       mcstBoruvka

  # Based on a cone-wise decomposition
    # Cone-wise rule:                     mcstCone



#### Bird rule ####

#' Bird's Rule for MCSTP
#'
#' @description
#' This function computes Bird's cost allocation (Bird, 1976) for a minimum cost spanning tree
#' problem \eqn{(N_0, C)}. The Bird rule, denoted as \eqn{B(N_0, C)}, assigns to each agent the
#' cost of the arc that connects them to the source during the construction of
#' a minimum cost spanning tree using Prim's algorithm.
#'
#' @details
#' For a given tree \eqn{g^n} constructed following Prim's algorithm, the allocation for
#' each agent \eqn{i \in N} is given by:
#'
#' \deqn{B_i(N_0, C) = c_{i^0 i},}
#'
#' where \eqn{i^0} is the node to which agent \eqn{i} is first connected in the process
#' of building the network.
#'
#' When the minimum cost spanning tree is not unique, the rule is extended by
#' averaging the allocations over all possible permutations \eqn{\pi \in \Pi_N}
#' (\code{method = "exact"}):
#'
#' \deqn{B(N_0, C) = \displaystyle{\dfrac{1}{n!} \sum_{\pi \in \Pi_N} B^{\pi}(N_0, C)},}
#'
#' where \eqn{B^{\pi}(N_0, C)} is the allocation obtained by applying Prim's
#' algorithm to \eqn{(N_0, C)} and solving indifferences by selecting the first
#' agent given by the permutation \eqn{\pi}, as proposed by Dutta and Kar (2004).
#'
#' Alternatively, \code{method = "montecarlo"} estimates the average allocation
#' through random sampling using \eqn{n_{\text{sim}}} permutations. This is highly useful
#' when the number of agents is large, as it significantly reduces computation
#' time while providing a very accurate approximation:
#'
#' \deqn{B(N_0, C) \approx \displaystyle{\dfrac{1}{n_{\text{sim}}} \sum_{k=1}^{n_{\text{sim}}} B^{\pi_k}(N_0, C)},}
#'
#' where \eqn{\pi_k} represents each of the randomly sampled permutations.
#' For networks with few agents, the exact method is preferred to avoid
#' potential overestimations.
#'
#' @param C cost matrix between nodes. Accepts multiple formats (see
#' "Supported formats for \code{C}" below).
#' @param draw logical; if \code{TRUE}, plots the network highlighting an optimal
#' tree in red, indicating the stage at which each arc is added (in brackets)
#' and the cost allocated to each agent (in parentheses below the node). For \eqn{n \le 3}
#' with ties and \code{method = "exact"}, it also displays all possible trees.
#' @param which.plot character string indicating which plots to display (only for \eqn{n \le 3}
#' with ties and \code{method = "exact"}). If \code{"details"}, displays the detailed breakdown of all trees and
#' allocations according to the agents' entry orders; if \code{"main"}, the average allocation
#' is plotted. Default is \code{c("main", "details")}.
#' @param titles logical; if \code{TRUE} (default), adds a main title
#' specifying the allocation rule and a subtitle with the algorithm, method (exact or Monte Carlo),
#' and the total network cost.
#' @param method character string specifying the calculation method when ties
#' occur. If \code{"exact"} (default), computes all \eqn{n!} permutations; if \code{"montecarlo"},
#' estimates the average allocation through random sampling using \code{nsim} permutations
#' (recommended for large networks).
#' @param nsim integer; the number of permutations to sample when using
#' \code{method = "montecarlo"}. Default is 10,000.
#'
#' @section Supported formats for \code{C}:
#'
#' Note: In all cases, \code{Inf} is accepted for disconnected node pairs. The first
#' row/column is assumed to be the source node (0).
#'
#' The function accepts the following formats for the cost input:
#'
#' \describe{
#'   \item{Numeric vector}{The lower triangle of a symmetric cost matrix, excluding the diagonal.
#'   Length must be \eqn{n(n+1)/2}, where \eqn{n} is the number of agents (excluding the source).}
#'
#'   \item{Square matrix / Data frame}{A standard adjacency matrix where \code{C[i, j]}
#'   represents the cost between node \eqn{i} and node \eqn{j}.}
#'
#'   \item{Edge list (matrix or data frame)}{An alternative to the full matrix.
#'   It must contain columns named \code{"from"} and \code{"to"},
#'   plus a third column for the weights/costs.}
#'
#'   \item{'igraph' object}{A weighted undirected graph. See the \code{\link[igraph:igraph-package]{igraph}}
#'   package for details on creating these objects.}
#'
#'   \item{File path}{A string pointing to a \code{.csv}, \code{.xlsx}, or \code{.xls} file.
#'   The file can contain either a square matrix (without headers) or an edge list
#'   (with \code{"from"} and \code{"to"} headers).}
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{bird}: the Bird's allocation vector \eqn{B(N_0, C)} for the natural ordering.
#'   \item \code{arcs}: a data frame of edges \eqn{(i, j)} in the final tree and the stage at which they were added.
#'   \item \code{e_bird}: the extended Bird's allocation (average over permutations).
#'   \item \code{total}: the total cost of the MCST, \eqn{m(N_0, C)}.
#'   \item \code{percentage}: the share of the total cost allocated to each agent.
#'   \item \code{ranking}: a ranking of agents by cost (from highest to lowest; ties marked with *).
#'   \item \code{nperms}: the number of permutations (\eqn{n!} or \code{nsim}) computed for the average allocation (if ties exist).
#'   \item \code{perms}: a data frame with the detailed allocations for the computed permutations (if ties exist).
#'   \item \code{is_unique}: logical; if \code{TRUE}, the MCST is unique.
#'   \item \code{method}: the method applied for computing the extended rule (\code{"exact"} or \code{"montecarlo"}).
#' }
#'
#' @seealso
#' \code{\link{mcstDuttaKar}} for other rule based on Prim's algorithm.
#'
#' \code{\link{mcstGameIrred}} for a cooperative game whose Shapley
#' value coincides with Bird's solution when applied to the irreducible matrix \eqn{\bar{C}}.
#'
#' \code{\link{mcstRules}} for an overview of the available rules and analysis tools in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' Bird C (1976) On cost allocation for a spanning tree: a game theoretic
#' approach. Networks 6(4):335–350.
#'
#' Dutta B, Kar A (2004) Cost monotonicity, consistency and minimum cost spanning tree
#' games. Games Econom Behav 48(2):223–248.
#'
#' @examples
#' # Simple vector input
#' mcstBird(c(12, 15, 20, 4, 6, 8), draw = TRUE)
#'
#' # Input with infinite costs (disconnected nodes)
#' C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
#'            Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
#' mcstBird(C_inf)
#'
#' # Matrix input with ties
#' C_mat <- matrix(c(0, 12, 15, 12,
#'                  12,  0,  4,  6,
#'                  15,  4,  0,  8,
#'                  12,  6,  8,  0), byrow = TRUE, ncol = 4)
#' mcstBird(C_mat, draw = TRUE, which.plot = "main", titles = FALSE)
#'
#' @concept Algorithmic Rules
#' @concept MCSTP
#'
#' @export

mcstBird <- function(C, draw = FALSE, which.plot = c("main", "details"), titles = TRUE, method = c("exact", "montecarlo"), nsim = 10000) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1      # Number of agents n = |N|
  N <- as.character(1:n)    # Set of agents N = {1, ..., n}

  method <- match.arg(method)
  which.plot <- match.arg(which.plot, several.ok = TRUE)

  # Inner function to calculate Prim's algorithm and Bird's allocation
  # for a given order pi
  .compute_bird <- function(C_mat, pi, n_agents) {
    S <- "0"          # Set of connected nodes. S^0 = {0} is the initial set (source)
    N_minus_S <- pi   # The set N \ S of agents not yet connected
    allocations <- numeric(n_agents)
    names(allocations) <- pi

    # Arcs (i, j) of the minimal spanning tree g^n
    arcs <- data.frame(
      i = character(n_agents),     # Origin node in S
      j = character(n_agents),     # New agent in N \ S
      stage = integer(n_agents),   # The stage p in which the arc was added
      stringsAsFactors = FALSE
    )

    tie_detected <- FALSE

    # Prim's Algorithm implementation
    for (p in 1:n_agents) {

      # Find arc (i, j) with i in S and j in N \ S with minimal cost c_ij
      sub_C <- C_mat[S, N_minus_S, drop = FALSE]
      min_cost <- min(sub_C)
      res <- which(sub_C == min_cost, arr.ind = TRUE)

      # Tie-breaking logic using the permutation pi
      if (nrow(res) > 1) {
        tie_detected <- TRUE

        # Select the first agent given by pi to solve indifferences
        candidates_j <- N_minus_S[res[, 2]]
        idx <- match(candidates_j, pi)
        best_idx <- which.min(idx)

        i <- S[res[best_idx, 1]]
        j <- N_minus_S[res[best_idx, 2]]

      } else {

        # Unique minimum cost arc (i, j)
        i <- S[res[1, 1]]
        j <- N_minus_S[res[1, 2]]
      }

      # Bird's allocation: agent j pays the connection cost to S
      # B_i = c_i*0_i
      allocations[j] <- min_cost
      arcs[p, ] <- list(i, j, p)

      # Update sets for the next stage p+1
      S <- c(S, j)
      N_minus_S <- setdiff(N_minus_S, j)
    }

    # Sort the allocation vector by agent index for the final output
    allocations <- allocations[order(as.numeric(names(allocations)))]
    return(list(allocations = allocations, arcs = arcs, tie_detected = tie_detected))
  }


  ## Execution ##

  # Standard Bird's Rule (using natural order of N)
  results <- .compute_bird(C_mat, N, n)
  m_cost <- sum(results$allocations) # Total cost m(N_0, C)

  # Extended Bird's Rule
  if (!results$tie_detected) {
    avg_allocation <- results$allocations
  } else {

    # Check if Monte Carlo is actually necessary
    if (method == "montecarlo" && factorial(n) <= nsim) {
      message("Warning: n! <= nsim. Consider using method = 'exact' for maximum efficiency")
    }

    # Performance validation
    if (method == "exact" && n > 10) {
      message("Warning: n > 10. Consider using method = 'montecarlo' to reduce processing time")
    }

    if (method == "exact") {
      # Generate all pi in Pi_N
      Pi_N <- gtools::permutations(n = n, r = n, v = N)
      nperms <- nrow(Pi_N)

      results_list <- lapply(1:nperms, function(k) {
        .compute_bird(C_mat, as.character(Pi_N[k, ]), n)
      })

      B_pi <- do.call(rbind, lapply(results_list, `[[`, "allocations"))
      arcs_pi <- lapply(results_list, `[[`, "arcs")
      perms_output <- Pi_N

    } else {

      nperms <- nsim

      results_list <- lapply(1:nsim, function(k) {
        current_pi <- sample(N)
        res <- .compute_bird(C_mat, current_pi, n)
        list(allocations = res$allocations, arcs = res$arcs, pi = current_pi)
      })

      B_pi <- do.call(rbind, lapply(results_list, `[[`, "allocations"))
      arcs_pi <- lapply(results_list, `[[`, "arcs")
      perms_output <- do.call(rbind, lapply(results_list, `[[`, "pi"))
    }

    # Calculate the average allocation
    avg_allocation <- colMeans(B_pi)

    perms <- data.frame(
      pi = apply(perms_output, 1, paste, collapse = "-"),
      B_pi,
      row.names = NULL,
      check.names = FALSE
    )
  }


  ## Visualization ##

  if (draw) {

    show_details <- "details" %in% which.plot
    show_main <- "main" %in% which.plot

    if (results$tie_detected && n <= 3 && method == "exact") {

      if (show_details) {
        rows <- if(n == 3) 2 else 1
        cols <- if(n == 3) 3 else 2
        old_par <- par(mfrow = c(rows, cols), oma = c(0, 0, if(titles) 5 else 0, 0))
        on.exit(par(old_par), add = TRUE)

        for (k in 1:nperms) {
          order <- perms$pi[k]
          current_allocation <- B_pi[k, ]
          .plot_mcstp(C_mat, arcs_pi[[k]], current_allocation,
                      main_title = "", sub_title = paste("Order:", order))
        }

        if (titles) {
          mtext("Detailed Bird Rule Allocations",
                side = 3, line = 1.5, cex = 1.2, font = 2, outer = TRUE)
          mtext(paste0("Algorithm: Prim  |  Method: exact  |  Total Cost: ", round(m_cost, 2)),
                side = 3, line = 0.1, cex = 0.7, family = "sans", font = 3,
                col = "#444444", outer = TRUE)
        }

        par(old_par)
      }

      if (show_details && show_main) {
        devAskNewPage(TRUE)
        on.exit(devAskNewPage(FALSE), add = TRUE)
      }

      if (show_main) {
        tit <- if(titles) "Extended Bird Rule Allocation" else ""
        sub_tit <- if(titles) paste0("Algorithm: Prim  |  Method: exact  |  Total Cost: ", round(m_cost, 2)) else ""
        .plot_mcstp(C_mat, results$arcs, avg_allocation,
                    main_title = tit, sub_title = sub_tit)
      }

    } else {
      tit <- if(titles) { if(results$tie_detected) "Extended Bird Rule Allocation" else "Bird Rule Allocation" } else ""
      lab <- if (method == "exact") {"exact"} else {"Monte Carlo"}
      sub_tit <- if(titles) paste0("Algorithm: Prim  |  Method: ", lab, "  |  Total Cost: ", round(m_cost, 2)) else ""
      final_alloc <- if(results$tie_detected) avg_allocation else results$allocations
      .plot_mcstp(C_mat, results$arcs, final_alloc,
                  main_title = tit, sub_title = sub_tit)
    }
  }

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$allocations, names(results$allocations))
  vals <- results$allocations[ord]; noms <- names(results$allocations)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    bird = results$allocations,
    arcs = results$arcs,
    e_bird = if(results$tie_detected) avg_allocation else results$allocations,
    total = m_cost,
    percentage = round((results$allocations / m_cost) * 100, 2),
    ranking = rank_star,
    nperms = if(results$tie_detected) nperms else NULL,
    perms = if(results$tie_detected) perms else NULL,
    is_unique = !results$tie_detected,
    method = method
  )

  class(output) <- "mcstp_bird"

  return(output)
}

#' @export
print.mcstp_bird <- function(x, ...) {
  if (x$is_unique) {
    print(round(x$bird, 2))
  } else {
    message("Non-unique MCST detected")
    alloc <- rbind(
      "B(id)" = x$bird,
      "E[B(pi)]" = x$e_bird
    )
    print(round(alloc, 2))
  }
  invisible(x)
}



#### Dutta-Kar rule ####

#' Dutta-Kar's Rule for MCSTP
#'
#' @description
#' This function computes Dutta-Kar's cost allocation (Dutta and Kar, 2004) for a minimum cost
#' spanning tree problem \eqn{(N_0, C)}. The Dutta-Kar rule, denoted as \eqn{DK(N_0, C)},
#' is defined through Prim's algorithm with a pivotal switch in the allocation
#' cost at each step to satisfy cost monotonicity.
#'
#' @details
#' Following Prim's algorithm, let \eqn{x^0 = 0}.
#' At each stage \eqn{p \in \{1, \dots, n\}}, the algorithm selects the minimum cost arc
#' \eqn{(i^p, j^p)} that connects a new agent \eqn{j^p} to the already connected set of
#' nodes. We update the accumulated cost as:
#'
#' \deqn{x^p = \max \{x^{p-1}, c_{i^p j^p}\}.}
#'
#' For each agent \eqn{i \in N}, there exists a stage \eqn{p(i)} such that
#' \eqn{i = j^{p(i)}}. Dutta-Kar's allocation for agent \eqn{i \in N} is defined as:
#'
#' \deqn{DK_i(N_0, C) = \min \{x^{p(i)-1}, c_{i^{p(i)} j^{p(i)}}\}.}
#'
#' When the minimum cost spanning tree is not unique, the rule is extended by
#' averaging the allocations over all possible permutations \eqn{\pi \in \Pi_N}
#' (\code{method = "exact"}):
#'
#' \deqn{DK(N_0, C) = \displaystyle{\frac{1}{n!} \sum_{\pi \in \Pi_N} DK^{\pi}(N_0, C)},}
#'
#' where \eqn{DK^{\pi}(N_0, C)} is the allocation obtained by applying Prim's algorithm
#' to \eqn{(N_0, C)} and solving indifferences by selecting the first agent given
#' by the permutation \eqn{\pi}.
#'
#' Alternatively, \code{method = "montecarlo"} estimates the average allocation
#' through random sampling using \eqn{n_{\text{sim}}} permutations. This is highly useful
#' when the number of agents is large, as it significantly reduces computation
#' time while providing a very accurate approximation:
#'
#' \deqn{DK(N_0, C) \approx \displaystyle{\dfrac{1}{n_{\text{sim}}} \sum_{k=1}^{n_{\text{sim}}} DK^{\pi_k}(N_0, C)},}
#'
#' where \eqn{\pi_k} represents each of the randomly sampled permutations.
#' For networks with few agents, the exact method is preferred to avoid
#' potential overestimations.
#'
#' @inheritParams mcstBird
#'
#' @inheritSection mcstBird Supported formats for \code{C}
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{dk}: the Dutta-Kar's allocation vector \eqn{DK(N_0, C)} for the natural ordering.
#'   \item \code{arcs}: a data frame of edges \eqn{(i, j)} in the final tree and the stage at which they were added.
#'   \item \code{e_dk}: the extended Dutta-Kar's allocation (average over permutations).
#'   \item \code{total}: the total cost of the MCST, \eqn{m(N_0, C)}.
#'   \item \code{percentage}: the share of the total cost allocated to each agent.
#'   \item \code{ranking}: a ranking of agents by cost (from highest to lowest; ties marked with *).
#'   \item \code{nperms}: the number of permutations (\eqn{n!} or \code{nsim}) computed for the average allocation (if ties exist).
#'   \item \code{perms}: a data frame with the detailed allocations for the computed permutations (if ties exist).
#'   \item \code{is_unique}: logical; if \code{TRUE}, the MCST is unique.
#'   \item \code{method}: the method applied for computing the extended rule (\code{"exact"} or \code{"montecarlo"}).
#' }
#'
#' @seealso
#' \code{\link{mcstBird}} for other rule based on Prim's algorithm.
#'
#' \code{\link{mcstRules}} for an overview of the available rules and analysis tools in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' Dutta B, Kar A (2004) Cost monotonicity, consistency and minimum cost spanning tree
#' games. Games Econom Behav 48(2):223–248
#'
#' @examples
#' # Simple vector input
#' mcstDuttaKar(c(12, 15, 20, 4, 6, 8), draw = TRUE)
#'
#' # Input with infinite costs (disconnected nodes)
#' C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
#'            Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
#' mcstDuttaKar(C_inf)
#'
#' # Matrix input with ties
#' C_mat <- matrix(c(0, 12, 15, 12,
#'                  12,  0,  4,  6,
#'                  15,  4,  0,  8,
#'                  12,  6,  8,  0), byrow = TRUE, ncol = 4)
#' mcstDuttaKar(C_mat, draw = TRUE, which.plot = "main", titles = FALSE)
#'
#' @concept Algorithmic Rules
#' @concept MCSTP
#'
#' @export

mcstDuttaKar <- function(C, draw = FALSE, which.plot = c("main", "details"), titles = TRUE, method = c("exact", "montecarlo"), nsim = 10000) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1      # Number of agents n = |N|
  N <- as.character(1:n)    # Set of agents N = {1, ..., n}

  method <- match.arg(method)
  which.plot <- match.arg(which.plot, several.ok = TRUE)

  # Inner function to calculate Prim's algorithm and Dutta-Kar's allocation
  # for a given order pi
  .compute_dk <- function(C_mat, pi, n_agents) {
    S <- "0"          # Set of connected nodes. S^0 = {0} is the initial set (source)
    N_minus_S <- pi   # The set N \ S of agents not yet connected
    allocations <- numeric(n_agents)
    names(allocations) <- pi

    # Arcs (i, j) of the minimal spanning tree g^n
    arcs <- data.frame(
      i = character(n_agents),     # Origin node in S
      j = character(n_agents),     # New agent in N \ S
      stage = integer(n_agents),   # The stage p in which the arc was added
      stringsAsFactors = FALSE
    )

    tie_detected <- FALSE

    # Trackers for Dutta-Kar
    j_p <- character(n_agents)     # Sequence of added agents j^p
    c_p <- numeric(n_agents)       # Connection costs c_{i^p j^p}
    x_p <- numeric(n_agents + 1)   # Values x^p
    x_p[1] <- 0                    # Stage 0: x^0 = 0

    # Prim's Algorithm implementation
    for (p in 1:n_agents) {

      # Find arc (i, j) with i in S and j in N \ S with minimal cost c_ij
      sub_C <- C_mat[S, N_minus_S, drop = FALSE]
      min_cost <- min(sub_C)
      res <- which(sub_C == min_cost, arr.ind = TRUE)

      # Tie-breaking logic using the permutation pi
      if (nrow(res) > 1) {
        tie_detected <- TRUE

        # Select the first agent given by pi to solve indifferences
        candidates_j <- N_minus_S[res[, 2]]
        idx <- match(candidates_j, pi)
        best_idx <- which.min(idx)

        i <- S[res[best_idx, 1]]
        j <- N_minus_S[res[best_idx, 2]]

      } else {

        # Unique minimum cost arc (i, j)
        i <- S[res[1, 1]]
        j <- N_minus_S[res[1, 2]]
      }

      # Dutta-Kar updates for stage p
      j_p[p] <- j
      c_p[p] <- min_cost
      x_p[p + 1] <- max(x_p[p], min_cost) # x^{p} = max{x^{p-1}, c_p}

      arcs[p, ] <- list(i, j, p)

      # Update sets for the next stage p+1
      S <- c(S, j)
      N_minus_S <- setdiff(N_minus_S, j)
    }

    # Dutta-Kar's allocation
    # DK_i = min{x^{p(i)-1}, c_{i^{p(i)} j^{p(i)}}}
    agent_stage <- setNames(1:n_agents, j_p)

    allocations <- sapply(pi, function(agent) {
      p_i <- agent_stage[agent]  # Stage in which agent was connected
      if (p_i < n_agents) {
        return(min(x_p[p_i + 1], c_p[p_i + 1]))
      } else {
        # Last agent pays the final accumulated value x^n
        return(x_p[n_agents + 1])
      }
    })

    # Sort the allocation vector by agent index for the final output
    allocations <- allocations[order(as.numeric(names(allocations)))]
    return(list(allocations = allocations, arcs = arcs, tie_detected = tie_detected))
  }


  ## Execution ##

  # Standard Dutta-Kar's Rule (using natural order of N)
  results <- .compute_dk(C_mat, N, n)
  m_cost <- sum(results$allocations) # Total cost m(N_0, C)

  # Extended Dutta-Kar's Rule
  if (!results$tie_detected) {
    avg_allocation <- results$allocations
  } else {

    # Check if Monte Carlo is actually necessary
    if (method == "montecarlo" && factorial(n) <= nsim) {
      message("Warning: n! <= nsim. Consider using method = 'exact' for maximum efficiency")
    }

    # Performance validation
    if (method == "exact" && n > 10) {
      message("Warning: n > 10. Consider using method = 'montecarlo' to reduce processing time")
    }

    if (method == "exact") {
      # Generate all pi in Pi_N
      Pi_N <- gtools::permutations(n = n, r = n, v = N)
      nperms <- nrow(Pi_N)

      results_list <- lapply(1:nperms, function(k) {
        .compute_dk(C_mat, as.character(Pi_N[k, ]), n)
      })

      DK_pi <- do.call(rbind, lapply(results_list, `[[`, "allocations"))
      arcs_pi <- lapply(results_list, `[[`, "arcs")
      perms_output <- Pi_N

    } else {

      nperms <- nsim

      results_list <- lapply(1:nsim, function(k) {
        current_pi <- sample(N)
        res <- .compute_dk(C_mat, current_pi, n)
        list(allocations = res$allocations, arcs = res$arcs, pi = current_pi)
      })

      DK_pi <- do.call(rbind, lapply(results_list, `[[`, "allocations"))
      arcs_pi <- lapply(results_list, `[[`, "arcs")
      perms_output <- do.call(rbind, lapply(results_list, `[[`, "pi"))
    }

    # Calculate the average allocation
    avg_allocation <- colMeans(DK_pi)

    perms <- data.frame(
      pi = apply(perms_output, 1, paste, collapse = "-"),
      DK_pi,
      row.names = NULL,
      check.names = FALSE
    )
  }


  ## Visualization ##

  if (draw) {

    show_details <- "details" %in% which.plot
    show_main <- "main" %in% which.plot

    if (results$tie_detected && n <= 3 && method == "exact") {

      if (show_details) {
        rows <- if(n == 3) 2 else 1
        cols <- if(n == 3) 3 else 2
        old_par <- par(mfrow = c(rows, cols), oma = c(0, 0, if(titles) 5 else 0, 0))
        on.exit(par(old_par), add = TRUE)

        for (k in 1:nperms) {
          order <- perms$pi[k]
          current_allocation <- DK_pi[k, ]
          .plot_mcstp(C_mat, arcs_pi[[k]], current_allocation,
                      main_title = "", sub_title = paste("Order:", order))
        }

        if (titles) {
          mtext("Detailed Dutta-Kar Rule Allocations",
                side = 3, line = 1.5, cex = 1.2, font = 2, outer = TRUE)
          mtext(paste0("Algorithm: Prim  |  Method: exact  |  Total Cost: ", round(m_cost, 2)),
                side = 3, line = 0.1, cex = 0.7, family = "sans", font = 3,
                col = "#444444", outer = TRUE)
        }

        par(old_par)
      }

      if (show_details && show_main) {
        devAskNewPage(TRUE)
        on.exit(devAskNewPage(FALSE), add = TRUE)
      }

      if (show_main) {
        tit <- if(titles) "Extended Dutta-Kar Rule Allocation" else ""
        sub_tit <- if(titles) paste0("Algorithm: Prim  |  Method: exact  |  Total Cost: ", round(m_cost, 2)) else ""
        .plot_mcstp(C_mat, results$arcs, avg_allocation,
                    main_title = tit, sub_title = sub_tit)
      }

    } else {
      tit <- if(titles) { if(results$tie_detected) "Extended Dutta-Kar Rule Allocation" else "Dutta-Kar Rule Allocation" } else ""
      lab <- if (method == "exact") {"exact"} else {"Monte Carlo"}
      sub_tit <- if(titles) paste0("Algorithm: Prim  |  Method: ", lab, "  |  Total Cost: ", round(m_cost, 2)) else ""
      final_alloc <- if(results$tie_detected) avg_allocation else results$allocations
      .plot_mcstp(C_mat, results$arcs, final_alloc,
                  main_title = tit, sub_title = sub_tit)
    }
  }

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$allocations, names(results$allocations))
  vals <- results$allocations[ord]; noms <- names(results$allocations)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    dk = results$allocations,
    arcs = results$arcs,
    e_dk = if(results$tie_detected) avg_allocation else results$allocations,
    total = m_cost,
    percentage = round((results$allocations / m_cost) * 100, 2),
    ranking = rank_star,
    nperms = if(results$tie_detected) nperms else NULL,
    perms = if(results$tie_detected) perms else NULL,
    is_unique = !results$tie_detected,
    method = method
  )

  class(output) <- "mcstp_dk"

  return(output)
}

#' @export
print.mcstp_dk <- function(x, ...) {
  if (x$is_unique) {
    print(round(x$dk, 2))
  } else {
    message("Non-unique MCST detected")
    alloc <- rbind(
      "DK(id)" = x$dk,
      "E[DK(pi)]" = x$e_dk
    )
    print(round(alloc, 2))
  }
  invisible(x)
}



#### Folk rule ####

#' Folk Rule for MCSTP
#'
#' @description
#' This function computes folk cost allocation (Feltkamp et al., 1994) for a
#' minimum cost spanning tree problem \eqn{(N_0, C)}. The folk rule, denoted as \eqn{F(N_0, C)}, also known
#' as the ERO (Equal Remaining Obligation) rule, is defined through
#' Kruskal's algorithm using a function with specific properties, called
#' the obligation function \eqn{o_i(S)}.
#'
#' @details
#' Following Kruskal's algorithm, let \eqn{g^p} be the network at stage \eqn{p \in \{1, \dots, n\}},
#' and \eqn{S_i^p := S(P(g^p), i)} the connected component containing agent \eqn{i}.
#' The folk rule is based on the obligation function:
#'
#' \deqn{o_i(S) = \dfrac{1}{|S|}.}
#'
#' At each stage \eqn{p}, let \eqn{c_{i^p j^p}} be the cost of the arc added to the
#' network. The folk allocation for each agent \eqn{i \in N} is given by:
#'
#' \deqn{F_i(N_0, C) = \displaystyle{ \sum_{p=1}^n c_{i^p j^p} \left( o_i(S_i^{p-1}) - o_i(S_i^p) \right)}.}
#'
#' @param C cost matrix between nodes. Accepts multiple formats (see
#' "Supported formats for \code{C}" below).
#' @param draw logical; if \code{TRUE}, plots the network highlighting an optimal
#' tree in red, indicating the stage at which each arc is added (in brackets)
#' and the cost allocated to each agent (in parentheses below the node).
#' @param titles logical; if \code{TRUE} (default), adds a main title
#' specifying the allocation rule and a subtitle with the algorithm used
#' and the total network cost.
#'
#' @inheritSection mcstBird Supported formats for \code{C}
#'
#' @note
#' The folk allocation is uniquely determined even when the minimum cost
#' spanning tree is not unique. Since the rule depends on the evolution of the
#' connected components rather than the specific arcs selected, any tie-breaking
#' selection in Kruskal's algorithm yields the same cost distribution.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{folk}: the folk allocation vector \eqn{F(N_0, C)}.
#'   \item \code{arcs}: a data frame of edges \eqn{(i, j)} in the final tree and the stage at which they were added.
#'   \item \code{total}: the total cost of the MCST, \eqn{m(N_0, C)}.
#'   \item \code{percentage}: the share of the total cost allocated to each agent.
#'   \item \code{ranking}: a ranking of agents by cost (from highest to lowest; ties marked with *).
#'   \item \code{is_unique}: logical; if \code{TRUE}, the MCST is unique.
#' }
#'
#' @seealso
#' \code{\link{mcstOWShapley}}, \code{\link{mcstPWShapley}} for other
#' rules based on Kruskal's algorithm.
#'
#' \code{\link{mcstBoruvka}} for an equivalent rule based on Boruvka's algorithm.
#'
#' \code{\link{mcstCone}} with \code{rule = "folk"} for the equivalent
#' implementation based on cone-wise decomposition.
#'
#' \code{\link{mcstGameIrred}}, \code{\link{mcstGameOpt}} for cooperative
#' games whose Shapley value coincides with the folk solution.
#'
#' \code{\link{mcstRules}} for an
#' overview of the available rules and analysis tools in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' Feltkamp V, Tijs S, Muto S (1994) On the irreducible core and the equal
#' remaining obligations rule of minimum cost spanning extension problems.
#' Technical Report 106, CentER DP 1994, Tilburg University, The Netherlands
#'
#' @examples
#' # Simple vector input
#' mcstFolk(c(12, 15, 20, 4, 6, 8), draw = TRUE)
#'
#' # Input with infinite costs (disconnected nodes)
#' C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
#'            Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
#' mcstFolk(C_inf)
#'
#' # Matrix input
#' C_mat <- matrix(c(0, 12, 15, 12,
#'                  12,  0,  4,  6,
#'                  15,  4,  0,  8,
#'                  12,  6,  8,  0), byrow = TRUE, ncol = 4)
#' mcstFolk(C_mat, draw = TRUE, titles = FALSE)
#'
#' @concept Algorithmic Rules
#' @concept MCSTP
#'
#' @export

mcstFolk <- function(C, draw = FALSE, titles = TRUE) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1      # Number of agents n = |N|
  N <- as.character(1:n)    # Set of agents N = {1, ..., n}

  # Inner function to calculate Kruskal's algorithm and folk allocation
  .compute_folk <- function(C_mat, n_agents, N_names) {
    N_0 <- c("0", N_names) # Set of nodes N_0 = N U {0}

    # Identify all possible arcs (i, j) and their costs c_ij
    # A^0(C) = {(i,j) | i,j in N_0, i != j}
    pairs <- combn(N_0, 2)
    A_0 <- data.frame(
      i = pairs[1, ],
      j = pairs[2, ],
      cost = apply(pairs, 2, function(x) C_mat[x[1], x[2]]),
      stringsAsFactors = FALSE
    )

    # Sort arcs by cost
    A_0 <- A_0[order(A_0$cost), ]

    tie_detected <- any(duplicated(A_0$cost))

    # Initial partition P: each node is in its own component S
    P <- lapply(N_0, function(x) c(x))
    names(P) <- N_0

    allocations <- numeric(n_agents)
    names(allocations) <- N_names

    # Arcs selected for the minimum spanning tree g^|N|
    g_N <- data.frame(i = character(n_agents), j = character(n_agents),
                      stage = integer(n_agents), stringsAsFactors = FALSE)

    p <- 1 # Stage counter

    # Kruskal's Algorithm implementation: sequential addition of arcs without introducing cycles
    for (m in 1:nrow(A_0)) {
      i <- A_0$i[m]
      j <- A_0$j[m]
      c_p <- A_0$cost[m] # Cost of the arc selected in stage p

      # Find current components S_i and S_j for nodes i and j
      idx_i <- which(sapply(P, function(x) i %in% x))
      idx_j <- which(sapply(P, function(x) j %in% x))

      # If nodes belong to different components, adding (i,j) doesn't form a cycle
      if (idx_i != idx_j) {
        S_i <- P[[idx_i]]
        S_j <- P[[idx_j]]
        source_in_i <- "0" %in% S_i
        source_in_j <- "0" %in% S_j


        # Folk allocation: o_i(S) = 1/|S|
        size_Si <- length(S_i)
        size_Sj <- length(S_j)
        size_total <- size_Si + size_Sj

        if (!source_in_i && !source_in_j) { # Both groups are not connected to the source
          allocations[S_i] <- allocations[S_i] + c_p * ((1/size_Si) - (1/size_total))
          allocations[S_j] <- allocations[S_j] + c_p * ((1/size_Sj) - (1/size_total))

        } else if (source_in_i && !source_in_j) { # S_i has the source
          allocations[S_j] <- allocations[S_j] + c_p * (1/size_Sj)

        } else if (!source_in_i && source_in_j) { # S_j has the source
          allocations[S_i] <- allocations[S_i] + c_p * (1/size_Si)
        }

        # Update sets for the next stage p+1
        g_N[p, ] <- list(i, j, p)
        p <- p + 1

        # Merge components to update partition P
        P[[idx_i]] <- c(S_i, S_j)
        P[[idx_j]] <- NULL
      }

      # Process is completed in |N| stages
      if (p > n_agents) break
    }

    # Sort the allocation vector by agent index for the final output
    allocations <- allocations[order(as.numeric(names(allocations)))]
    return(list(allocations = allocations, arcs = g_N, tie_detected = tie_detected))
  }


  ## Execution ##

  results <- .compute_folk(C_mat, n, N)
  m_cost <- sum(results$allocations) # Total cost m(N_0, C)


  ## Visualization ##

  if (draw) {

    tit <- if(titles) "Folk Rule Allocation" else ""
    sub_tit <- if(titles) paste0("Algorithm: Kruskal  |  Total Cost: ", round(m_cost, 2)) else ""
    .plot_mcstp(C_mat, results$arcs, results$allocations,
                main_title = tit, sub_title = sub_tit)
  }

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$allocations, names(results$allocations))
  vals <- results$allocations[ord]; noms <- names(results$allocations)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    folk = results$allocations,
    arcs = results$arcs,
    total = m_cost,
    percentage = round((results$allocations / m_cost) * 100, 2),
    ranking = rank_star,
    is_unique = !results$tie_detected
  )

  class(output) <- "mcstp_folk"

  return(output)
}

#' @export
print.mcstp_folk <- function(x, ...) {
  print(round(x$folk, 2))
  invisible(x)
}



#### Optimistic weighted Shapley rule ####

#' Optimistic Weighted Shapley's Rule for MCSTP
#'
#' @description
#' This function computes the optimistic weighted Shapley cost allocation (Bergantiños and Lorenzo- Freire, 2008a,b) for a
#' minimum cost spanning tree problem \eqn{(N_0, C)}. The optimistic weighted Shapley rule,
#' denoted as \eqn{f^{\varrho^{ow}}(N_0, C)}, is defined through Kruskal's algorithm
#' using a function with specific properties, called the obligation function \eqn{o_i(S)}.
#'
#' @details
#' Following Kruskal's algorithm, let \eqn{g^p} be the network at stage \eqn{p \in \{1, \dots, n\}},
#' and \eqn{S_i^p := S(P(g^p), i)} the connected component containing agent \eqn{i}.
#' The optimistic weighted Shapley rule is based on the obligation function:
#'
#' \deqn{o_i(S) = \dfrac{w_i}{\sum_{j \in S} w_j},}
#'
#' where \eqn{w_i > 0} is the weight associated with agent \eqn{i}.
#'
#' At each stage \eqn{p}, let \eqn{c_{i^p j^p}} be the cost of the arc added to the
#' network. The optimistic weighted Shapley allocation for each agent \eqn{i \in N} is given by:
#'
#' \deqn{f^{\varrho^{ow}}_i(N_0, C) = \displaystyle{ \sum_{p=1}^n c_{i^p j^p} \left( o_i(S_i^{p-1}) - o_i(S_i^p) \right)}.}
#'
#' @inheritParams mcstFolk
#' @param weights a numeric vector of strictly positive weights for each agent in \eqn{N}.
#' The length of the vector must match the number of agents.
#'
#' @inheritSection mcstBird Supported formats for \code{C}
#'
#' @note
#' The optimistic weighted Shapley allocation is uniquely determined even when the minimum
#' cost spanning tree is not unique. Since the rule depends on the evolution of the
#' connected components rather than the specific arcs selected, any tie-breaking
#' selection in Kruskal's algorithm yields the same cost distribution.
#'
#' Furthermore, when all agents share the same weight, the optimistic weighted Shapley
#' allocation is mathematically equivalent to the \code{\link{mcstFolk}}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{ows}: the optimistic weighted Shapley allocation vector \eqn{f^{\varrho^{ow}}(N_0, C)}.
#'   \item \code{weights}: the weight vector used for the allocation.
#'   \item \code{arcs}: a data frame of edges \eqn{(i, j)} in the final tree and the stage at which they were added.
#'   \item \code{total}: the total cost of the MCST, \eqn{m(N_0, C)}.
#'   \item \code{percentage}: the share of the total cost allocated to each agent.
#'   \item \code{ranking}: a ranking of agents by cost (from highest to lowest; ties marked with *).
#'   \item \code{is_unique}: logical; if \code{TRUE}, the MCST is unique.
#' }
#'
#' @seealso
#' \code{\link{mcstFolk}}, \code{\link{mcstPWShapley}} for other
#' rules based on Kruskal's algorithm.
#'
#' \code{\link{mcstCone}} with \code{rule = "owshapley"} for the equivalent
#' implementation based on cone-wise decomposition.
#'
#' \code{\link{mcstRules}} for an
#' overview of the available rules and analysis tools in the package.
#'
#' @references
#' Bergantiños G, Lorenzo-Freire S (2008a) A characterization of optimistic
#' weighted Shapley rules in mini- mum cost spanning tree problems. Econ
#' Theor 35(3):523–538
#'
#' Bergantiños G, Lorenzo-Freire S (2008b) “Optimistic” weighted Shapley
#' rules in minimum cost spanning tree problems. European J Operatl Res
#' 185(1):289–298
#'
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' @examples
#' # Simple vector input with custom weights
#' mcstOWShapley(c(12, 15, 20, 4, 6, 8), weights = c(2, 1, 1))
#'
#' # Matrix input with equal weights (equivalent to folk rule)
#' C_mat <- matrix(c(0, 12, 15, 12,
#'                  12,  0,  4,  6,
#'                  15,  4,  0,  8,
#'                  12,  6,  8,  0), byrow = TRUE, ncol = 4)
#' mcstOWShapley(C_mat, weights = rep(1,3), draw = TRUE, titles = FALSE)
#'
#' @concept Algorithmic Rules
#' @concept MCSTP
#'
#' @export

mcstOWShapley <- function(C, weights, draw = FALSE, titles = TRUE) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1      # Number of agents n = |N|
  N <- as.character(1:n)    # Set of agents N = {1, ..., n}

  # Weight validation
  if (is.null(names(weights))) {
    if (length(weights) != n) {
      stop("Weight vector length must match the number of agents")
    }
    names(weights) <- N
  }

  # Inner function to calculate Kruskal's algorithm and Optimistic weighted Shapley's allocation
  .compute_ows <- function(C_mat, w, n_agents, N_names) {
    N_0 <- c("0", N_names) # Set of nodes N_0 = N U {0}

    # Identify all possible arcs (i, j) and their costs c_ij
    # A^0(C) = {(i,j) | i,j in N_0, i != j}
    pairs <- combn(N_0, 2)
    A_0 <- data.frame(
      i = pairs[1, ],
      j = pairs[2, ],
      cost = apply(pairs, 2, function(x) C_mat[x[1], x[2]]),
      stringsAsFactors = FALSE
    )

    # Sort arcs by cost
    A_0 <- A_0[order(A_0$cost), ]

    tie_detected <- any(duplicated(A_0$cost))

    # Initial partition P: each node is in its own component S
    P <- lapply(N_0, function(x) c(x))
    names(P) <- N_0

    allocations <- numeric(n_agents)
    names(allocations) <- N_names

    # Arcs selected for the minimum spanning tree g^|N|
    g_N <- data.frame(i = character(n_agents), j = character(n_agents),
                      stage = integer(n_agents), stringsAsFactors = FALSE)

    p <- 1 # Stage counter

    # Kruskal's Algorithm implementation: sequential addition of arcs without introducing cycles
    for (m in 1:nrow(A_0)) {
      i <- A_0$i[m]
      j <- A_0$j[m]
      c_p <- A_0$cost[m] # Cost of the arc selected in stage p

      # Find current components S_i and S_j for nodes i and j
      idx_i <- which(sapply(P, function(x) i %in% x))
      idx_j <- which(sapply(P, function(x) j %in% x))

      # If nodes belong to different components, adding (i,j) doesn't form a cycle
      if (idx_i != idx_j) {
        S_i <- P[[idx_i]]
        S_j <- P[[idx_j]]
        source_in_i <- "0" %in% S_i
        source_in_j <- "0" %in% S_j


        # Optimistic Weighted Shapley's allocation: o_i(S) = w_i / sum(w_j for j in S)
        S_i_agents <- S_i[S_i != "0"]
        S_j_agents <- S_j[S_j != "0"]

        W_i <- sum(w[S_i_agents])
        W_j <- sum(w[S_j_agents])
        W_total <- W_i + W_j

        if (!source_in_i && !source_in_j) { # Both groups are not connected to the source
          allocations[S_i_agents] <- allocations[S_i_agents] + c_p * w[S_i_agents] * ((1/W_i) - (1/W_total))
          allocations[S_j_agents] <- allocations[S_j_agents] + c_p * w[S_j_agents] * ((1/W_j) - (1/W_total))

        } else if (source_in_i && !source_in_j) { # S_i has the source
          allocations[S_j_agents] <- allocations[S_j_agents] + c_p * (w[S_j_agents] / W_j)

        } else if (!source_in_i && source_in_j) { # S_j has the source
          allocations[S_i_agents] <- allocations[S_i_agents] + c_p * (w[S_i_agents] / W_i)
        }

        # Update sets for the next stage p+1
        g_N[p, ] <- list(i, j, p)
        p <- p + 1

        # Merge components to update partition P
        P[[idx_i]] <- c(S_i, S_j)
        P[[idx_j]] <- NULL
      }

      # Process is completed in |N| stages
      if (p > n_agents) break
    }

    # Sort the allocation vector by agent index for the final output
    allocations <- allocations[order(as.numeric(names(allocations)))]
    return(list(allocations = allocations, arcs = g_N, tie_detected = tie_detected))
  }


  ## Execution ##

  results <- .compute_ows(C_mat, weights, n, N)
  m_cost <- sum(results$allocations) # Total cost m(N_0, C)


  ## Visualization ##

  if (draw) {

    tit <- if(titles) "Optimistic Weighted Shapley Rule Allocation" else ""
    sub_tit <- if(titles) paste0("Algorithm: Kruskal  |  Total Cost: ", round(m_cost, 2)) else ""
    .plot_mcstp(C_mat, results$arcs, results$allocations,
                main_title = tit, sub_title = sub_tit)
  }

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$allocations, names(results$allocations))
  vals <- results$allocations[ord]; noms <- names(results$allocations)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    ows = results$allocations,
    weights = weights,
    arcs = results$arcs,
    total = m_cost,
    percentage = round((results$allocations / m_cost) * 100, 2),
    ranking = rank_star,
    is_unique = !results$tie_detected
  )

  class(output) <- "mcstp_ows"

  return(output)
}

#' @export
print.mcstp_ows <- function(x, ...) {
  cat(paste0("Weights: ", paste(x$weights, collapse = " "), "\n"))
  print(round(x$ows, 2))
  invisible(x)
}



#### Pessimistic weighted Shapley rule ####

#' Pessimistic Weighted Shapley's Rule for MCSTP
#'
#' @description
#' This function computes the pessimistic weighted Shapley cost allocation (Lorenzo and Lorenzo-Freire, 2009) for a
#' minimum cost spanning tree problem \eqn{(N_0, C)}. The pessimistic weighted Shapley rule,
#' denoted as \eqn{f^{\varrho^{pw}}(N_0, C)}, is defined through Kruskal's algorithm
#' using a function with specific properties, called the obligation function \eqn{o_i(S)}.
#'
#' @details
#' Following Kruskal's algorithm, let \eqn{g^p} be the network at stage \eqn{p \in \{1, \dots, n\}},
#' and \eqn{S_i^p := S(P(g^p), i)} the connected component containing agent \eqn{i}.
#' The pessimistic weighted Shapley rule is based on the obligation function:
#'
#' \deqn{o_i(S) = \displaystyle{ \sum_{\pi \in \Pi(S \setminus \{i\})} \prod_{j=1}^{|S|-1} \dfrac{w_{\pi^{-1}(j)}}{\sum_{k=1}^{j} w_{\pi^{-1}(k)} + w_i},}}
#'
#' where \eqn{\Pi(\cdot)} denotes the set of all permutations
#' of a given set, and \eqn{w_i > 0} is the weight associated with agent \eqn{i}.
#'
#' At each stage \eqn{p}, let \eqn{c_{i^p j^p}} be the cost of the arc added to the
#' network. The pessimistic weighted Shapley allocation for each agent \eqn{i \in N} is given by:
#'
#' \deqn{f^{\varrho^{pw}}_i(N_0, C) = \displaystyle{ \sum_{p=1}^n c_{i^p j^p} \left( o_i(S_i^{p-1}) - o_i(S_i^p) \right)}.}
#'
#' @inheritParams mcstOWShapley
#'
#' @inheritSection mcstBird Supported formats for \code{C}
#'
#' @note
#' The pessimistic weighted Shapley allocation is uniquely determined even when the minimum
#' cost spanning tree is not unique. Since the rule depends on the evolution of the
#' connected components rather than the specific arcs selected, any tie-breaking
#' selection in Kruskal's algorithm yields the same cost distribution.
#'
#' Furthermore, when all agents share the same weight, the pessimistic weighted Shapley
#' allocation is mathematically equivalent to the \code{\link{mcstFolk}}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{pws}: the pessimistic weighted Shapley allocation vector \eqn{f^{\varrho^{pw}}(N_0, C)}.
#'   \item \code{weights}: the weight vector used for the allocation.
#'   \item \code{arcs}: a data frame of edges \eqn{(i, j)} in the final tree and the stage at which they were added.
#'   \item \code{total}: the total cost of the MCST, \eqn{m(N_0, C)}.
#'   \item \code{percentage}: the share of the total cost allocated to each agent.
#'   \item \code{ranking}: a ranking of agents by cost (from highest to lowest; ties marked with *).
#'   \item \code{is_unique}: logical; if \code{TRUE}, the MCST is unique.
#' }
#'
#' @seealso
#' \code{\link{mcstFolk}}, \code{\link{mcstOWShapley}} for other
#' rules based on Kruskal's algorithm.
#'
#' \code{\link{mcstRules}} for an
#' overview of the available rules and analysis tools in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' Lorenzo L, Lorenzo-Freire S (2009) A characterization of Kruskal sharing rules
#' for minimum cost spanning tree problems. Int J Game Theory 38(1):107–126
#'
#' @examples
#' # Simple vector input with custom weights
#' mcstPWShapley(c(12, 15, 20, 4, 6, 8), weights = c(2, 1, 1))
#'
#' # Matrix input with equal weights (equivalent to folk rule)
#' C_mat <- matrix(c(0, 12, 15, 12,
#'                  12,  0,  4,  6,
#'                  15,  4,  0,  8,
#'                  12,  6,  8,  0), byrow = TRUE, ncol = 4)
#' mcstPWShapley(C_mat, weights = rep(1,3), draw = TRUE, titles = FALSE)
#'
#' @concept Algorithmic Rules
#' @concept MCSTP
#'
#' @export

mcstPWShapley <- function(C, weights, draw = FALSE, titles = TRUE) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1      # Number of agents n = |N|
  N <- as.character(1:n)    # Set of agents N = {1, ..., n}

  # Weight validation
  if (is.null(names(weights))) {
    if (length(weights) != n) {
      stop("Weight vector length must match the number of agents")
    }
    names(weights) <- N
  }

  # Inner function to calculate Kruskal's algorithm and Pessimistic weighted Shapley's allocation
  .compute_pws <- function(C_mat, w, n_agents, N_names) {
    N_0 <- c("0", N_names) # Set of nodes N_0 = N U {0}

    # Helper function for pessimistic obligation o_i(S)
    .compute_o_pess <- function(i, S, w_vec) {
      if (length(S) == 1 && S[1] == i) return(1)
      S_minus_i <- setdiff(S, i)

      f <- function(rem, current_sum) {
        if (length(rem) == 0) return(1)
        res <- 0
        for (k in seq_along(rem)) {
          next_w <- w_vec[rem[k]]
          new_sum <- current_sum + next_w
          term <- next_w / new_sum
          res <- res + term * f(rem[-k], new_sum)
        }
        return(res)
      }
      return(f(S_minus_i, w_vec[i]))
    }

    # Identify all possible arcs (i, j) and their costs c_ij
    # A^0(C) = {(i,j) | i,j in N_0, i != j}
    pairs <- combn(N_0, 2)
    A_0 <- data.frame(
      i = pairs[1, ],
      j = pairs[2, ],
      cost = apply(pairs, 2, function(x) C_mat[x[1], x[2]]),
      stringsAsFactors = FALSE
    )

    # Sort arcs by cost
    A_0 <- A_0[order(A_0$cost), ]

    tie_detected <- any(duplicated(A_0$cost))

    # Initial partition P: each node is in its own component S
    P <- lapply(N_0, function(x) c(x))
    names(P) <- N_0

    allocations <- numeric(n_agents)
    names(allocations) <- N_names

    # Arcs selected for the minimum spanning tree g^|N|
    g_N <- data.frame(i = character(n_agents), j = character(n_agents),
                      stage = integer(n_agents), stringsAsFactors = FALSE)

    p <- 1 # Stage counter

    # Kruskal's Algorithm implementation: sequential addition of arcs without introducing cycles
    for (m in 1:nrow(A_0)) {
      i <- A_0$i[m]
      j <- A_0$j[m]
      c_p <- A_0$cost[m] # Cost of the arc selected in stage p

      # Find current components S_i and S_j for nodes i and j
      idx_i <- which(sapply(P, function(x) i %in% x))
      idx_j <- which(sapply(P, function(x) j %in% x))

      # If nodes belong to different components, adding (i,j) doesn't form a cycle
      if (idx_i != idx_j) {
        S_i <- P[[idx_i]]
        S_j <- P[[idx_j]]
        source_in_i <- "0" %in% S_i
        source_in_j <- "0" %in% S_j

        # Pessimistic Weighted Shapley's allocation: based on product-sum obligation o_i(S)
        agents_i <- S_i[S_i != "0"]
        agents_j <- S_j[S_j != "0"]
        agents_union <- c(agents_i, agents_j)

        if (!source_in_i && !source_in_j) { # Both groups are not connected to the source
          o_old_i <- sapply(agents_i, .compute_o_pess, S = agents_i, w_vec = w)
          o_new_i <- sapply(agents_i, .compute_o_pess, S = agents_union, w_vec = w)
          allocations[agents_i] <- allocations[agents_i] + c_p * (o_old_i - o_new_i)

          o_old_j <- sapply(agents_j, .compute_o_pess, S = agents_j, w_vec = w)
          o_new_j <- sapply(agents_j, .compute_o_pess, S = agents_union, w_vec = w)
          allocations[agents_j] <- allocations[agents_j] + c_p * (o_old_j - o_new_j)

        } else if (source_in_i && !source_in_j) { # S_i has the source
          o_old_j <- sapply(agents_j, .compute_o_pess, S = agents_j, w_vec = w)
          allocations[agents_j] <- allocations[agents_j] + c_p * o_old_j

        } else if (!source_in_i && source_in_j) { # S_j has the source
          o_old_i <- sapply(agents_i, .compute_o_pess, S = agents_i, w_vec = w)
          allocations[agents_i] <- allocations[agents_i] + c_p * o_old_i
        }

        # Update sets for the next stage p+1
        g_N[p, ] <- list(i, j, p)
        p <- p + 1

        # Merge components to update partition P
        P[[idx_i]] <- c(S_i, S_j)
        P[[idx_j]] <- NULL
      }

      # Process is completed in |N| stages
      if (p > n_agents) break
    }

    # Sort the allocation vector by agent index for the final output
    allocations <- allocations[order(as.numeric(names(allocations)))]
    return(list(allocations = allocations, arcs = g_N, tie_detected = tie_detected))
  }


  ## Execution ##

  results <- .compute_pws(C_mat, weights, n, N)
  m_cost <- sum(results$allocations) # Total cost m(N_0, C)


  ## Visualization ##

  if (draw) {

    tit <- if(titles) "Pessimistic Weighted Shapley Rule Allocation" else ""
    sub_tit <- if(titles) paste0("Algorithm: Kruskal  |  Total Cost: ", round(m_cost, 2)) else ""
    .plot_mcstp(C_mat, results$arcs, results$allocations,
                main_title = tit, sub_title = sub_tit)
  }

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$allocations, names(results$allocations))
  vals <- results$allocations[ord]; noms <- names(results$allocations)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    pws = results$allocations,
    weights = weights,
    arcs = results$arcs,
    total = m_cost,
    percentage = round((results$allocations / m_cost) * 100, 2),
    ranking = rank_star,
    is_unique = !results$tie_detected
  )

  class(output) <- "mcstp_pws"

  return(output)
}

#' @export
print.mcstp_pws <- function(x, ...) {
  cat(paste0("Weights: ", paste(x$weights, collapse = " "), "\n"))
  print(round(x$pws, 2))
  invisible(x)
}


#### Boruvka rule ####

#' Boruvka's Rule for MCSTP
#'
#' @description
#' This function computes Boruvka's cost allocation (Bergantiños and Vidal- Puga, 2011) for a minimum
#' cost spanning tree problem \eqn{(N_0, C)}. Boruvka's rule, denoted as
#' \eqn{\beta^{\pi}(N_0, C)}, is defined through Boruvka's algorithm.
#'
#' @details
#' The Boruvka allocation divides the cost of the optimal tree through a step-by-step
#' process over \eqn{\gamma} stages. At each stage \eqn{s}, agents pay the maximum
#' possible proportion \eqn{p^s} of the cheapest arc associated with their connected
#' component.
#'
#' Let \eqn{A^s} be the set of non-completely paid arcs, \eqn{\varrho_{ij}^s} the
#' proportion of the cost of arc \eqn{(i,j)} already paid, and \eqn{N_{ij}^s} the set of agents
#' paying for arc \eqn{(i,j)} in stage \eqn{s}. The proportion \eqn{p^s} is given by:
#'
#' \deqn{p^s = \min \left\{ \dfrac{1 - \varrho_{ij}^{s-1}}{|N_{ij}^s|} : (i,j) \in A^{s-1}, ~ N_{ij}^s \neq \emptyset \right\}.}
#'
#' Let \eqn{a_i^s} be the arc partially paid by agent \eqn{i} in stage \eqn{s}. The cost paid
#' by agent \eqn{i} in this stage is:
#'
#' \deqn{f_i^s = p^s c_{a_i^s}.}
#'
#' The process finishes in \eqn{\gamma} stages when the tree is completely paid, i.e.,
#' \eqn{\sum_{s=1}^\gamma p^s = 1}. The Boruvka allocation for each agent
#' \eqn{i \in N} is given by:
#'
#' \deqn{\beta_i^{\pi}(N_0, C) = \displaystyle{\sum_{s=1}^\gamma f_i^s.}}
#'
#' @inheritParams mcstFolk
#'
#' @inheritSection mcstBird Supported formats for \code{C}
#'
#' @note
#' Boruvka's allocation is uniquely determined even when the minimum cost
#' spanning tree is not unique. Since the rule depends on the evolution of the
#' connected components rather than the specific arcs selected, any tie-breaking
#' selection in Boruvka's algorithm yields the same cost distribution.
#'
#' Furthermore, as demonstrated by Bergantiños and Vidal-Puga (2011), for any
#' order \eqn{\pi} used in the algorithm, Boruvka's allocation is
#' equivalent to the \code{\link{mcstFolk}}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{boruvka}: the Boruvka's allocation vector \eqn{\beta^{\pi}(N_0, C)}.
#'   \item \code{arcs}: a data frame of edges \eqn{(i, j)} in the final tree and the stage at which they were added.
#'   \item \code{total}: the total cost of the MCST, \eqn{m(N_0, C)}.
#'   \item \code{percentage}: the share of the total cost allocated to each agent.
#'   \item \code{ranking}: a ranking of agents by cost (from highest to lowest; ties marked with *).
#'   \item \code{is_unique}: logical; if \code{TRUE}, the MCST is unique.
#' }
#'
#' @seealso
#' \code{\link{mcstFolk}} for an equivalent rule based on Kruskal's algorithm.
#'
#' \code{\link{mcstRules}} for an overview of the available rules and analysis tools in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2011) The folk solution and Boruvka’s algorithm
#' in minimum cost spanning tree problems. Discrete Appl Math 159(12):1279–1283
#'
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' @examples
#' # Simple vector input
#' mcstBoruvka(c(12, 15, 20, 4, 6, 8), draw = TRUE)
#'
#' # Input with infinite costs (disconnected nodes)
#' C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
#'            Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
#' mcstBoruvka(C_inf)
#'
#' # Matrix input
#' C_mat <- matrix(c(0, 12, 15, 12,
#'                  12,  0,  4,  6,
#'                  15,  4,  0,  8,
#'                  12,  6,  8,  0), byrow = TRUE, ncol = 4)
#' mcstBoruvka(C_mat, draw = TRUE, titles = FALSE)
#'
#' @concept Algorithmic Rules
#' @concept MCSTP
#'
#' @export

mcstBoruvka <- function(C, draw = FALSE, titles = TRUE) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1      # Number of agents n = |N|
  N <- as.character(1:n)    # Set of agents N = {1, ..., n}

  # Inner function to calculate Boruvka's algorithm and allocation
  .compute_boruvka <- function(C_mat, n_agents, N_names) {
    N_0 <- c("0", N_names) # Set of nodes N_0 = N U {0}

    # Identify all possible arcs (i, j) and their costs c_ij
    pairs <- combn(N_0, 2)
    A_0 <- data.frame(
      i = pairs[1, ],
      j = pairs[2, ],
      cost = apply(pairs, 2, function(x) C_mat[x[1], x[2]]),
      stringsAsFactors = FALSE
    )

    # Sort arcs by cost
    A_0 <- A_0[order(A_0$cost, A_0$i, A_0$j), ]

    tie_detected <- any(duplicated(A_0$cost))

    # Arcs selected for the minimum spanning tree g^|N|
    g_N <- data.frame(i = character(), j = character(),
                      cost = numeric(), stringsAsFactors = FALSE)

    P_map <- setNames(N_0, N_0)
    find_P <- function(node) {
      while(P_map[node] != node)
        node <- P_map[node]
      node
    }

    # Boruvka's Algorithm implementation: merging components until only one remains
    while(length(unique(sapply(N_0, find_P))) > 1) {
      m_s <- c() # Arcs selected in this stage s
      roots <- unique(sapply(N_0, find_P))

      for (r in roots) {
        if (r == find_P("0")) next
        comp_nodes <- N_0[sapply(N_0, find_P) == r]

        # Cheapest arc connecting T and N_0 \ T
        idx <- which(((A_0$i %in% comp_nodes & !(A_0$j %in% comp_nodes)) |
                        (A_0$j %in% comp_nodes & !(A_0$i %in% comp_nodes))))

        if (length(idx) > 0)
          m_s <- c(m_s, idx[1])
      }

      for (idx in unique(m_s)) {
        i <- A_0$i[idx]
        j <- A_0$j[idx]

        if (find_P(i) != find_P(j)) {
          g_N <- rbind(g_N, A_0[idx, ])
          P_map[find_P(j)] <- find_P(i)
        }
      }
    }

    # Re-order tree arcs to follow order pi for the allocation rule
    g_N <- g_N[order(match(paste(g_N$i, g_N$j), paste(A_0$i, A_0$j))), ]

    # Boruvka's allocation
    e_ij <- setNames(rep(0, nrow(g_N)), 1:nrow(g_N)) # Paid proportion of each arc
    allocations <- setNames(rep(0, n), N)

    while(any(e_ij < 0.9999)) {
      # Partition based on completely paid arcs (e_ij = 1)
      P_paid <- setNames(N_0, N_0)
      find_paid <- function(node) {
        while(P_paid[node] != node)
          node <- P_paid[node]
        node
      }
      idx_paid <- which(e_ij >= 0.9999)

      if(length(idx_paid) > 0) {
        for(idx in idx_paid) {
          r_i <- find_paid(g_N$i[idx])
          r_j <- find_paid(g_N$j[idx])
          if(r_i != r_j)
            P_paid[r_j] <- r_i
        }
      }

      # Identify agents paying for each arc in this stage (N_ij)
      N_ij <- list()
      roots <- unique(sapply(N_0, find_paid))
      root_0 <- find_paid("0")

      for (r in roots) {
        if (r == root_0) next
        comp_agents <- N[sapply(N, find_paid) == r]

        for(idx in 1:nrow(g_N)) {
          if ((find_paid(g_N$i[idx]) == r && find_paid(g_N$j[idx]) != r) ||
              (find_paid(g_N$j[idx]) == r && find_paid(g_N$i[idx]) != r)) {
            N_ij[[as.character(idx)]] <- unique(c(N_ij[[as.character(idx)]], comp_agents))
            break
          }
        }
      }

      # Maximum proportion p_s that can be paid in this stage
      p_s <- 1.0
      for (k in names(N_ij)) {
        val <- (1 - e_ij[k]) / length(N_ij[[k]])

        if (val < p_s)
          p_s <- val
      }

      # Update allocations and payment proportions
      for (k in names(N_ij)) {
        idx <- as.numeric(k)
        for (agent in N_ij[[k]]) allocations[agent] <- allocations[agent] + p_s * g_N$cost[idx]
        e_ij[k] <- e_ij[k] + length(N_ij[[k]]) * p_s
      }
    }

    tree_arcs <- data.frame(i = g_N$i, j = g_N$j, stage = 1:nrow(g_N))

    # Sort the allocation vector by agent index for the final output
    allocations <- allocations[order(as.numeric(names(allocations)))]
    return(list(allocations = allocations, arcs = tree_arcs, tie_detected = tie_detected))
  }


  ## Execution ##

  results <- .compute_boruvka(C_mat, n, N)
  m_cost <- sum(results$allocations) # Total cost m(N_0, C)


  ## Visualization ##

  if (draw) {

    tit <- if(titles) "Boruvka Rule Allocation" else ""
    sub_tit <- if(titles) paste0("Algorithm: Boruvka  |  Total Cost: ", round(m_cost, 2)) else ""
    .plot_mcstp(C_mat, results$arcs, results$allocations,
                main_title = tit, sub_title = sub_tit)
  }

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$allocations, names(results$allocations))
  vals <- results$allocations[ord]; noms <- names(results$allocations)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    boruvka = results$allocations,
    arcs = results$arcs,
    total = m_cost,
    percentage = round((results$allocations / m_cost) * 100, 2),
    ranking = rank_star,
    is_unique = !results$tie_detected
  )

  class(output) <- "mcstp_boruvka"

  return(output)
}

#' @export
print.mcstp_boruvka <- function(x, ...) {
  print(round(x$boruvka, 2))
  invisible(x)
}



#### Cone-wise decomposition ####

#' Cone-wise Decomposition Rules for MCSTP
#'
#' @description
#' This function computes cost allocations for minimum cost spanning tree problems
#' \eqn{(N_0, C)} using the cone-wise decomposition introduced by Norde et al. (2004). This
#' technique allows extending rules defined for elementary MCSTPs to any
#' MCSTP.
#'
#' @details
#' Any MCSTP can be written as a nonnegative combination of elementary
#' MCST \eqn{C = \sum_{q=1}^{m(C)} x^q C^q}, where the costs of the arcs in
#' \eqn{C^q} are 0 or 1. A rule \eqn{R} is extended as:
#'
#' \deqn{R(N_0, C) = \displaystyle{\sum_{q=1}^{m(C)} x^q R(N_0, C^q).}}
#'
#' Let \eqn{g^q} be the network at stage \eqn{q \in \{1, \dots, m(C)\}}, and
#' \eqn{S_i^q := S(P(g^q), i)} the connected component containing agent \eqn{i} in
#' the graph induced by zero-cost arcs.
#' For each elementary problem \eqn{(N_0, C^q)}, the agents in a component \eqn{S}
#' not containing the source (0) divide the cost of connecting to the source
#' according to the chosen rule:
#' \itemize{
#'   \item Folk rule (\code{"folk"}): \deqn{F(N_0, C^q) = \begin{cases}
#'    \dfrac{1}{|S_i^q|} & \text{if } 0 \notin S_i^q \\
#'    0 & \text{otherwise}.
#'    \end{cases}}
#'
#'   \item Optimistic weighted Shapley rule (\code{"owshapley"}): \deqn{f^{\varrho^{ow}}(N_0, C^q) = \begin{cases}
#'    \dfrac{w_i}{\sum_{j \in S_i^q} w_j} & \text{if } 0 \notin S_i^q \\
#'    0 & \text{otherwise}.
#'    \end{cases}}
#'
#'   \item Bogomolnaia and Moulin family (\code{"bogomolnaia"}):
#'
#'   For \eqn{\lambda \in [0, +\infty)},
#'    \deqn{R^\lambda(N_0, C^q) = \begin{cases}
#'    \dfrac{\lambda^{\delta_i}}{\sum_{j \in S_i^q} \lambda^{\delta_j}} & \text{if } 0 \notin S_i^q \\
#'    0 & \text{otherwise},
#'    \end{cases}}
#'   where \eqn{\delta_i} denotes the number of non-null arcs in \eqn{C^q} containing
#'   agent \eqn{i}.
#'
#'   For \eqn{\lambda = +\infty},
#'    \deqn{R^\lambda(N_0, C^q) = \begin{cases}
#'    \arg\max_{j \in S_i^q} \delta_j & \text{if } 0 \notin S_i^q \\
#'    0 & \text{otherwise},
#'    \end{cases}}
#'    i.e., agents attaining the maximum value of \eqn{\delta_j} share the cost equally.
#' }
#'
#' @inheritParams mcstFolk
#' @param rule a character string indicating the rule \eqn{R} to be applied to the elementary MCSTPs.
#' One of \code{"folk"} (default), \code{"owshapley"}, or \code{"bogomolnaia"}.
#' @param weights a numeric vector of strictly positive weights for each agent in \eqn{N}.
#' The length of the vector must match the number of agents; only if \code{rule = "owshapley"}.
#' @param lambda a non-negative numeric parameter (\code{1} by default); only if \code{rule = "bogomolnaia"}.
#'
#' @inheritSection mcstBird Supported formats for \code{C}
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{allocations}: the cost allocation vector \eqn{R(N_0, C)} for the chosen rule \eqn{R}.
#'   \item \code{rule}: the rule \eqn{R} applied for the allocation.
#'   \item \code{total}: the total cost of the MCST, \eqn{m(N_0, C)}.
#'   \item \code{percentage}: the share of the total cost allocated to each agent.
#'   \item \code{ranking}: a ranking of agents by cost (from highest to lowest; ties marked with *).
#'   \item \code{decomposition}: a matrix showing the allocation at each cost level of the decomposition.
#'   \item \code{weights}: the weight vector used for the allocation; only if \code{rule = "owshapley"}.
#'   \item \code{lambda}: the lambda value used for the allocation; only if \code{rule = "bogomolnaia"}.
#' }
#'
#' @note
#' This function focuses on the analytical and numerical decomposition of the cost
#' allocation. For graphical representations, use the direct rule functions
#' (e.g., \code{\link{mcstFolk}}).
#'
#' The allocation for \code{rule = "folk"} is also obtained when all weights are
#' equal in \code{rule = "owshapley"}, or when \eqn{\lambda = 1} in
#' \code{rule = "bogomolnaia"}.
#'
#' @seealso
#' \code{\link{mcstFolk}}, \code{\link{mcstOWShapley}} for direct implementations of these rules.
#'
#' \code{\link{mcstRules}} for an overview of the available rules and analysis tools in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' Norde H, Moretti S, Tijs S (2004) Minimum cost spanning tree games and population
#' monotonic allocation schemes. Eur J Oper Res 154(1):84–97
#'
#' @examples
#' C <- c(12, 15, 20, 4, 6, 8)
#'
#' # Folk rule
#' mcstCone(C, rule = "folk")
#'
#' # Optimistic weighted Shapley rule
#' mcstCone(C, rule = "owshapley", weights = c(2, 1, 1))
#'
#' # Bogomolnaia and Moulin family
#' mcstCone(C, rule = "bogomolnaia", lambda = 2)
#' mcstCone(C, rule = "bogomolnaia", lambda = Inf)
#'
#' @concept Algorithmic Rules
#' @concept MCSTP
#'
#' @export

mcstCone <- function(C, rule = c("folk", "owshapley", "bogomolnaia"), weights = NULL, lambda = 1) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1      # Number of agents n = |N|
  N <- as.character(1:n)    # Set of agents N = {1, ..., n}

  rule <- match.arg(rule)

  # Rule and parameters validation
  if (rule == "owshapley") {
    if (is.null(weights)) stop("Weights are required for the 'owshapley' rule")

    if (length(weights) != n) stop("Weight vector length must match the number of agents")

    if (is.null(names(weights))) {
      names(weights) <- N
    }
  }

  if (rule == "bogomolnaia") {
    if (is.null(lambda) || is.na(lambda) || lambda < 0) {
      stop("Lambda must be a non-negative number")
    }
  }

  # Inner function to calculate cone-wise decomposition and allocation
  # R(N_0, C) = sum_{q=1}^{m(C)} x^q * R(N_0, C^q)
  .compute_conewise <- function(C_mat, rule_type, w, l, n_agents, N_names) {
    N_0 <- c("0", N) # Set of nodes N_0 = N U {0}

    # Identify unique costs to define the decomposition
    unique_costs <- sort(unique(as.vector(C_mat[upper.tri(C_mat)])))
    unique_costs <- unique_costs[unique_costs > 0 & is.finite(unique_costs)]

    # x^q: non-negative weights for the decomposition (differences between unique costs)
    x_q <- diff(c(0, unique_costs))

    allocations <- setNames(numeric(n_agents), N_names)
    step_details <- matrix(0, nrow = length(unique_costs), ncol = n_agents, dimnames = list(paste0("Cost_", unique_costs), N_names))

    # Iterate through each elementary mcstp (C^q)
    for (q in seq_along(unique_costs)) {

      # C^q: cost matrix where c_ij = 1 if c_ij >= unique_costs[q] and 0 otherwise
      C_q <- (C_mat >= unique_costs[q]) * 1
      diag(C_q) <- 0

      # For each elementary problem, we find components S connected by paths of cost 0
      adj_zero <- (C_q == 0)
      visited <- setNames(rep(FALSE, length(N_0)), N_0)
      elem_alloc <- setNames(numeric(n_agents), N_names)

      for (node in N_0) {
        if (!visited[node]) {
          # Component S Identification
          comp <- node
          stack <- node
          visited[node] <- TRUE

          while(length(stack) > 0) {
            curr <- stack[1]
            stack <- stack[-1]
            neighbors <- N_0[adj_zero[curr, ] & !visited]
            visited[neighbors] <- TRUE
            comp <- c(comp, neighbors)
            stack <- c(stack, neighbors)
          }

          # Agents connected to the source pay 0
          if (!("0" %in% comp)) {
            S <- comp[comp != "0"] # Set of agents in the component S

            # Application of the rule R to the elementary problem (N_0, C^q)
            if (rule_type == "folk") {
              # Folk rule: 1/|S|
              elem_alloc[S] <- 1 / length(S)

            } else if (rule_type == "owshapley") {
              # Optimistic weighted Shapley: w_i / sum(w_j for j in S)
              elem_alloc[S] <- w[S] / sum(w[S])

            } else if (rule_type == "bogomolnaia") {
              # Number of non-null arcs in C^q containing agent i (delta_i)
              deltas <- sapply(S, function(i) sum(C_q[i, ] == 1))
              if (is.infinite(l)) {
                # Lambda = infinity: split among agents with max delta
                is_max <- (deltas == max(deltas))
                elem_alloc[S[is_max]] <- 1 / sum(is_max)
              } else {
                # General lambda formula
                scores <- l^deltas
                if (sum(scores) == 0) { # Handling lambda = 0
                  is_min <- (deltas == min(deltas))
                  elem_alloc[S[is_min]] <- 1 / sum(is_min)
                } else {
                  elem_alloc[S] <- scores / sum(scores)
                }
              }
            }
          }
        }
      }

      # Accumulate results for the general mcstp
      step_res <- x_q[q] * elem_alloc
      step_details[q, ] <- step_res
      allocations <- allocations + step_res
    }

    # Sort the allocation vector by agent index for the final output
    allocations <- allocations[order(as.numeric(names(allocations)))]
    return(list(allocations = allocations, steps = step_details))
  }


  ## Execution ##

  results <- .compute_conewise(C_mat, rule, weights, lambda, n, N)
  m_cost <- sum(results$allocations) # Total cost m(N_0, C)

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$allocations, names(results$allocations))
  vals <- results$allocations[ord]; noms <- names(results$allocations)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    allocations = results$allocations,
    rule = rule,
    total = m_cost,
    percentage = round((results$allocations / m_cost) * 100, 2),
    ranking = rank_star,
    decomposition = results$steps
  )

  if (rule == "owshapley") output$weights <- weights
  if (rule == "bogomolnaia") output$lambda <- lambda

  class(output) <- "mcstp_conewise"

  return(output)
}

#' @export
print.mcstp_conewise <- function(x, ...) {
  if (x$rule == "owshapley") cat(paste0("Weights: ", paste(x$weights, collapse = " "), "\n"))
  if (x$rule == "bogomolnaia") cat(paste0("Lambda: ", x$lambda, "\n"))
  print(round(x$allocations, 2))
  invisible(x)
}



# Modified 12/05 21:00
