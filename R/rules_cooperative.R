# RULES DEFINED THROUGH COOPERATIVE GAMES

  # Private game:           private_game
  # Irreducible game:       irred_game
  # Optimistic game:        opt_game
  # Public game:            public_game
  # Cycle-complete game:    cc_game



#### Private Game ####

#' Cost Allocations based on the Private Game for MCSTP (Kar's Rule)
#'
#' @description
#' This function computes the cost allocation associated with the private game,
#' also known as the pessimistic game (Bird, 1976), for a minimum cost spanning
#' tree problem \eqn{(N_0, C)}. The private game, denoted as \eqn{(N, v_C^p)} \eqn{(}or simply
#' \eqn{(N, v^p)}\eqn{)}, evaluates the cost of each coalition \eqn{S \subseteq N} under the assumption
#' that agents outside the coalition are not available to provide their nodes
#' as intermediaries. When the Shapley value is applied as the solution concept to this game,
#' the resulting allocation is known as the Kar's rule (Kar, 2002).
#'
#' @details
#' The private game \eqn{(N, v_C^p)} is a cooperative game with
#' transferable utility (TU game) defined for each coalition
#' \eqn{S \subseteq N} as:
#'
#' \deqn{v_C^p(S) = m(S_0, C),}
#'
#' where \eqn{m(S_0, C)} is the cost of the minimum spanning tree of the
#' subproblem induced by the agents in \eqn{S} and the source.
#' This approach is termed "private" or "pessimistic" because the nodes in
#' \eqn{N \setminus S} belong to these agents, and their participation
#' is needed in order to use their nodes.
#'
#' Kar's rule is defined by applying the Shapley value to this TU
#' game, which assigns to each agent \eqn{i \in N} the average of
#' their marginal contributions over all possible permutations \eqn{\pi \in \Pi_N}:
#'
#' \deqn{Sh_i(N, v^p) = \displaystyle{\dfrac{1}{n!} \sum_{\pi \in \Pi_N} [v^p(Pre(i, \pi) \cup \{i\}) - v^p(Pre(i, \pi))]},}
#'
#' where \eqn{Pre(i, \pi)} is the set of players that precede agent \eqn{i}
#' in the ordering \eqn{\pi}.
#'
#' @param C a symmetric square matrix or a numeric vector representing the
#' lower triangle of costs (ordered by columns) among the nodes in \eqn{N_0}.
#' The first row and column are assumed to be the source (0). Supports \code{Inf}
#' for disconnected nodes.
#' @param sol character; the solution concept to apply to the game.
#' Currently, only \code{sol = "shapley"} (default) is supported.
#' @param draw logical; if \code{TRUE}, plots the network highlighting an optimal
#' tree in red (constructed using Prim's algorithm), indicating the cost allocated to
#' each agent in brackets below the node. For \eqn{n \le 3}, it also displays the
#' optimal network for each subcoalition \eqn{S \subset N}.
#' @param which numeric vector indicating which plots to display. If \code{1} (only
#' for \eqn{n \le 3}), displays the optimal network for each subcoalition
#' \eqn{S \subset N}; if \code{2}, the allocation for \eqn{N} is plotted.
#' Default is \code{c(1, 2)}.
#' @param titles logical; if \code{TRUE} (default), adds a main title
#' specifying the cooperative game and a subtitle with the solution concept used
#' and the total network cost.
#'
#' @note
#' The private game is a strong reference for stability.
#' Classical rules such as Bird, Dutta-Kar, or the folk rule are known to belong
#' to the core of \eqn{v^p} for any MCST problem. As a result, these
#' allocations are coalitionally stable under the pessimistic cost perspective.
#'
#' Beyond the Shapley value, other solution concepts have been studied
#' for this game, including the core and the nucleolus. It has been proven
#' that computing both the Shapley value and the nucleolus is NP-hard.
#'
#' @return A list containing:
#' \itemize{
#'    \item \code{coalitions}: a character vector representing all possible coalitions \eqn{S \subseteq N} ordered by cardinality.
#'    \item \code{v_p}: a vector containing the characteristic function values \eqn{v^p(S)} for all coalitions.
#'    \item \code{contributions}: a data frame with the marginal contributions for all permutations \eqn{\pi \in \Pi_N}.
#'    \item \code{shapley}: the Shapley value allocation vector (Kar's rule).
#'    \item \code{total}: the total cost of the MCST, \eqn{v^p(N)}.
#'    \item \code{nperms}: the number of permutations \eqn{n!} used to compute marginal contributions and the resulting Shapley value.
#'    \item \code{percentage}: the share of the total cost allocated to each agent.
#'    \item \code{ranking}: ranking of agents by cost (from highest to lowest; ties marked with *).
#'    }
#'
#' @seealso
#' \code{\link{bird_rule}}, \code{\link{dk_rule}}, \code{\link{folk_rule}}
#' for classical rules belonging to the core of \eqn{v^p}.
#'
#' \code{\link{irred_game}} for the pessimistic game based on the irreducible cost matrix \eqn{\bar{C}}.
#'
#' \code{\link{cc_game}} for the pessimistic game based on the cycle-complete cost matrix \eqn{C^*}.
#'
#' \code{\link{opt_game}}, \code{\link{public_game}} for optimistic alternatives to the private approach.
#'
#' \code{\link{alloc_rules}} for an overview of the available rules in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' Bird C (1976) On cost allocation for a spanning tree: a game theoretic
#' approach. Networks 6(4):335–350
#'
#' Kar A (2002) Axiomatization of the Shapley value on minimum cost
#' spanning tree games. Games Econom Behav 38(2):265–277
#'
#' @examples
#' # Simple vector input
#' private_game(c(12, 15, 20, 4, 6, 8), sol = "shapley", draw = TRUE)
#'
#' # Input with infinite costs (disconnected nodes)
#' C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
#'            Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
#' private_game(C_inf, sol = "shapley")
#'
#' # Matrix input
#' C_mat <- matrix(c(0, 12, 15, 12,
#'                  12,  0,  4,  6,
#'                  15,  4,  0,  8,
#'                  12,  6,  8,  0), byrow = TRUE, ncol = 4)
#' private_game(C_mat, draw = TRUE, which = 2, titles = FALSE)
#'
#' @concept Cooperative Games
#' @concept MCSTP
#'
#' @export

private_game <- function(C, sol = "shapley", draw = FALSE, which = c(1, 2), titles = TRUE) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1  # Number of agents n = |N|

  # Build the characteristic function v^p(S) using binary indexing
  ncoals <- 2^n - 1
  vvals <- numeric(ncoals)

  for (i in 1:ncoals) {
    # Decode binary index i to find members of coalition S
    S <- which(as.logical(intToBits(i)[1:n]))
    sub_idx <- c(1, S + 1) # Include the source (0)
    sub_C <- C_mat[sub_idx, sub_idx, drop = FALSE]

    # Compute MST cost for this specific subset
    vvals[i] <- .get_cost(sub_C)
  }


  ## Execution ##

  # Calculate the allocation based on the chosen solution concept
  if (sol == "shapley") {
    results <- .compute_shapley(vvals, n)
  } else {
    stop("Solution not implemented. Currently, only 'shapley' is available")
  }

  m_cost <- vvals[ncoals] # Total cost v^p(N)

  # Generate all possible coalitions ordered by cardinality for final output
  coalitions <- unlist(lapply(1:n, function(k) combn(n, k, simplify = FALSE)), recursive = FALSE)
  coals_names <- sapply(coalitions, function(S) paste0("{", paste(S, collapse = ","), "}"))
  coals_names[length(coals_names)] <- "N"
  idx <- sapply(coalitions, function(S) sum(2^(S - 1)))

  # To display characteristic function v^p(S)
  values <- vvals[idx]
  # v_S <- as.data.frame(t(values))
  # colnames(v_S) <- coals_names; rownames(v_S) <- "v^p(S)"
  # v_S <- data.frame(
  #   S = coals_names,
  #   "v^p(S)" = values,
  #   row.names = NULL,
  #   check.names = FALSE
  # )

  # cat("Private Game\n")
  # cat("Characteristic function v^p(S):\n")
  # print(v_S)
  # cat("---------------\n")
  # cat(paste0("Total Cost v^p(N): ", round(m_cost, 2), "\n"))
  # cat("\nSolution: Shapley value\n")
  # print(round(results$value, 2))

  # if (n <= 3) {
  #   cat("\nDisplaying marginal contributions:\n")
  #   nperms <- nrow(results$contributions)
  #   cat(paste0("Computed ", nperms, " permutations\n"))
  #   print(results$contributions)
  # } else {
  #   nperms <- factorial(n)
  #   cat(paste0("\nMarginal contributions: computed ", nperms, " permutations (too many to display individually)\n"))
  #   }


  ## Visualization ##

  if (draw) {

    if (is.null(which)) which <- c(1, 2)
    show_details <- 1 %in% which
    show_main <- 2 %in% which

    if (n > 1 && n <= 3 && show_details) {

      rows <- if(n == 3) 2 else 1
      cols <- if(n == 3) 3 else 2
      old_par <- par(mfrow = c(rows, cols), oma = c(0, 0, if(titles) 5 else 0, 0))
      on.exit(par(old_par), add = TRUE)

      for (i in 1:(ncoals - 1)) {
        S <- coalitions[[i]]
        sub_idx <- c(1, S + 1)

        C_inf <- matrix(Inf, nrow = n + 1, ncol = n + 1, dimnames = list(0:n, 0:n))

        C_inf[sub_idx, sub_idx] <- C_mat[sub_idx, sub_idx]
        sub_C <- C_mat[sub_idx, sub_idx, drop = FALSE]
        arcs_S <- .get_arcs(sub_C)

        .plot_mcstp(C_inf, arcs_S, NULL,
                    main_title = "", sub_title = paste0("S = {", paste(S, collapse=","), "}")) #   |  v^p(S): ", round(values[i], 2)
      }

      if (titles) {
        mtext("MCST for each subcoalition S",
              side = 3, line = 1.5, cex = 1.2, font = 2, outer = TRUE)
        mtext("Private Game (N, v^p)", side = 3, line = 0.1, cex = 0.7, family = "sans", font = 3,
              col = "#444444", outer = TRUE) # paste0("Private Game  |  Total Network Cost v^p(N): ", round(m_cost, 2))
      }

      par(old_par)
    }

    if (n > 1 && n <= 3 && show_details && show_main) {
      devAskNewPage(TRUE)
      on.exit(devAskNewPage(FALSE), add = TRUE)
    }

    if (show_main) {
      arcs_N <- .get_arcs(C_mat)
      tit <- if(titles) "Private Game (N, v^p)" else ""
      sub_tit <- if(titles) paste0("Solution: Shapley value (Kar's Rule)  |  Total Network Cost v^p(N): ", round(m_cost, 2)) else ""
      .plot_mcstp(C_mat, arcs_N, results$value,
                  main_title = tit, sub_title = sub_tit)
    }
  }

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$value, names(results$value))
  vals <- results$value[ord]; noms <- names(results$value)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    coalitions = coals_names,
    v_p = values,
    contributions = results$contributions,
    shapley = results$value,
    total = m_cost,
    nperms = nrow(results$contributions),
    percentage = round((results$value / m_cost) * 100, 2),
    ranking = rank_star
  )

  class(output) <- "mcstp_private"

  return(output)
}

#' @export
print.mcstp_private <- function(x, ...) {
  v_S <- data.frame(
    S = x$coalitions,
    "v^p(S)" = x$v_p,
    row.names = NULL,
    check.names = FALSE
  )
  print(v_S)
  cat("\nSolution: Shapley value\n")
  print(round(x$shapley, 2))
  invisible(x)
}



#### Irreducible Game ####

#' Cost Allocations based on the Irreducible Game for MCSTP
#'
#' @description
#' This function computes the cost allocation associated with the irreducible game
#' (Bird, 1976) for a minimum cost spanning tree problem \eqn{(N_0, C)}. The
#' irreducible game, denoted as \eqn{(N, v_C^i)} \eqn{(}or simply
#' \eqn{(N, v^i)}\eqn{)}, is defined as the private game associated with
#' the irreducible cost matrix \eqn{\bar{C}}. When the Shapley value is applied
#' as the solution concept to this game, the resulting allocation coincides
#' with the folk rule of \eqn{(N_0, C)}, and Bird's rule of \eqn{(N_0, \bar{C})}
#' (Bergantiños and Vidal-Puga, 2007a).
#'
#' @details
#' The irreducible game \eqn{(N, v_C^i)} is a cooperative game with
#' transferable utility (TU game) defined for each coalition
#' \eqn{S \subseteq N} as:
#'
#' \deqn{v_C^i(S) = m(S_0, \bar{C}),}
#'
#' where \eqn{m(S_0, \bar{C})} is the cost of the minimum spanning tree of the
#' subproblem induced by the agents in \eqn{S} and the source, evaluated over
#' the irreducible costs \eqn{\bar{C}}. The irreducible cost between any two nodes
#' corresponds to the maximum cost arc on the unique path connecting them in
#' the optimal tree.
#'
#' Applying the Shapley value to this TU game assigns to each agent \eqn{i \in N}
#' the average of their marginal contributions over all possible permutations
#' \eqn{\pi \in \Pi_N}:
#'
#' \deqn{Sh_i(N, v^i) = \displaystyle{\frac{1}{n!} \sum_{\pi \in \Pi_N} [v^i(Pre(i, \pi) \cup \{i\}) - v^i(Pre(i, \pi))]},}
#'
#' where \eqn{Pre(i, \pi)} is the set of players that precede agent \eqn{i}
#' in the ordering \eqn{\pi}.
#'
#' @inheritParams private_game
#'
#' @note
#' The irreducible game is a strong reference for stability. Its core (the
#' irreducible core) is always non-empty and constitutes a subset of the core
#' of \eqn{v^p}. Classical rules such as Bird and folk are
#' known to belong to the core of \eqn{v^i}.
#'
#' Unlike the private game, calculating solution concepts for the irreducible
#' game is computationally efficient. It has been proven that computing both
#' the Shapley value and the nucleolus of \eqn{v^i} can be done in \eqn{O(|N|^2)}
#' time.
#'
#' @return A list containing:
#' \itemize{
#'    \item \code{coalitions}: a character vector representing all possible coalitions \eqn{S \subseteq N} ordered by cardinality.
#'    \item \code{C_irred}: the irreducible cost matrix \eqn{\bar{C}} associated with the MCST of the problem.
#'    \item \code{v_i}: a vector containing the characteristic function values \eqn{v^i(S)} for all coalitions.
#'    \item \code{contributions}: a data frame with the marginal contributions for all permutations \eqn{\pi \in \Pi_N}.
#'    \item \code{shapley}: the Shapley value allocation vector.
#'    \item \code{total}: the total cost of the MCST, \eqn{v^i(N)}.
#'    \item \code{nperms}: the number of permutations \eqn{n!} used to compute marginal contributions and the resulting Shapley value.
#'    \item \code{percentage}: the share of the total cost allocated to each agent.
#'    \item \code{ranking}: ranking of agents by cost (from highest to lowest; ties marked with *).
#' }
#'
#' @seealso
#' \code{\link{bird_rule}} of \eqn{(N_0, \bar{C})}, \code{\link{folk_rule}}
#' of \eqn{(N_0, C)} for classical rules belonging to the core of \eqn{v^i}
#' that coincide with \eqn{Sh(N,v^i)}.
#'
#' \code{\link{private_game}}, \code{\link{public_game}} for the pessimistic and optimistic games based on the original costs.
#'
#' \code{\link{cc_game}} for the pessimistic game based on the cycle-complete cost matrix \eqn{C^*}.
#'
#' \code{\link{opt_game}} for an equivalent game under the Shapley value.
#'
#' \code{\link{alloc_rules}} for an overview of the available rules in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2007a) A fair rule in minimum cost spanning
#' tree problems. J Econom Theory 137(1):326–352
#'
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' Bird C (1976) On cost allocation for a spanning tree: a game theoretic
#' approach. Networks 6(4):335–350
#'
#' @examples
#' # Simple vector input
#' irred_game(c(12, 15, 20, 4, 6, 8), sol = "shapley", draw = TRUE)
#'
#' # Input with infinite costs (disconnected nodes)
#' C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
#'            Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
#' irred_game(C_inf, sol = "shapley")
#'
#' # Matrix input
#' C_mat <- matrix(c(0, 12, 15, 12,
#'                  12,  0,  4,  6,
#'                  15,  4,  0,  8,
#'                  12,  6,  8,  0), byrow = TRUE, ncol = 4)
#' irred_game(C_mat, draw = TRUE, which = 2, titles = FALSE)
#'
#' @concept Cooperative Games
#' @concept MCSTP
#'
#' @export

irred_game <- function(C, sol = "shapley", draw = FALSE, which = c(1, 2), titles = TRUE) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1  # Number of agents n = |N|

  # Compute the irreducible matrix barC
  C_irred <- .get_irreducible(C_mat)

  # Build the characteristic function v^i(S) using binary indexing over barC
  ncoals <- 2^n - 1
  vvals <- numeric(ncoals)

  for (i in 1:ncoals) {
    # Decode binary index i to find members of coalition S
    S <- which(as.logical(intToBits(i)[1:n]))
    sub_idx <- c(1, S + 1) # Include the source (node 0)
    sub_C <- C_irred[sub_idx, sub_idx, drop = FALSE]

    # Compute MST cost for this specific subset
    vvals[i] <- .get_cost(sub_C)
  }


  ## Execution ##

  # Calculate the allocation based on the chosen solution concept
  if (sol == "shapley") {
    results <- .compute_shapley(vvals, n)
  } else {
    stop("Solution not implemented. Currently, only 'shapley' is available")
  }

  m_cost <- vvals[ncoals] # Total cost v^i(N)

  # Generate all possible coalitions ordered by cardinality for final output
  coalitions <- unlist(lapply(1:n, function(k) combn(n, k, simplify = FALSE)), recursive = FALSE)
  coals_names <- sapply(coalitions, function(S) paste0("{", paste(S, collapse = ","), "}"))
  coals_names[length(coals_names)] <- "N"
  idx <- sapply(coalitions, function(S) sum(2^(S - 1)))

  # To display characteristic function v^i(S)
  values <- vvals[idx]
  # v_S <- as.data.frame(t(values))
  # colnames(v_S) <- coals_names; rownames(v_S) <- "v^i(S)"
  # v_S <- data.frame(
  #   S = coals_names,
  #   "v^i(S)" = values,
  #   row.names = NULL,
  #   check.names = FALSE
  # )

  # cat("Irreducible Game\n")
  # cat("Characteristic function v^i(S):\n")
  # print(v_S)
  # cat("---------------\n")
  # cat(paste0("Total Cost v^i(N): ", round(m_cost, 2), "\n"))
  # cat("\nSolution: Shapley value\n")
  # print(round(results$value, 2))

  # if (n <= 3) {
  #   cat("\nDisplaying marginal contributions:\n")
  #   nperms <- nrow(results$contributions)
  #   cat(paste0("Computed ", nperms, " permutations\n"))
  #   print(results$contributions)
  # } else {
  #   nperms <- factorial(n)
  #   cat(paste0("\nMarginal contributions: computed ", nperms, " permutations (too many to display individually)\n"))
  # }


  ## Visualization ##

  if (draw) {

    if (is.null(which)) which <- c(1, 2)
    show_details <- 1 %in% which
    show_main <- 2 %in% which

    if (n > 1 && n <= 3 && show_details) {

      rows <- if(n == 3) 2 else 1
      cols <- if(n == 3) 3 else 2
      old_par <- par(mfrow = c(rows, cols), oma = c(0, 0, if(titles) 5 else 0, 0))
      on.exit(par(old_par), add = TRUE)

      for (i in 1:(ncoals - 1)) {
        S <- coalitions[[i]]
        sub_idx <- c(1, S + 1)

        C_inf <- matrix(Inf, nrow = n + 1, ncol = n + 1, dimnames = list(0:n, 0:n))

        C_inf[sub_idx, sub_idx] <- C_irred[sub_idx, sub_idx]
        sub_C <- C_irred[sub_idx, sub_idx, drop = FALSE]
        arcs_S <- .get_arcs(sub_C)

        .plot_mcstp(C_inf, arcs_S, NULL,
                    main_title = "", sub_title = paste0("S = {", paste(S, collapse=","), "}")) #   |  v^i(S): ", round(values[i], 2)
      }

      if (titles) {
        mtext("MCST for each subcoalition S",
              side = 3, line = 1.5, cex = 1.2, font = 2, outer = TRUE)
        mtext("Irreducible Game (N, v^i)", side = 3, line = 0.1, cex = 0.7, family = "sans", font = 3,
              col = "#444444", outer = TRUE) # paste0("Irreducible Game  |  Total Network Cost v^i(N): ", round(m_cost, 2))
      }

      par(old_par)
    }

    if (n > 1 && n <= 3 && show_details && show_main) {
      devAskNewPage(TRUE)
      on.exit(devAskNewPage(FALSE), add = TRUE)
    }

    if (show_main) {
      arcs_N <- .get_arcs(C_irred)
      tit <- if(titles) "Irreducible Game (N, v^i)" else ""
      sub_tit <- if(titles) paste0("Solution: Shapley value  |  Total Network Cost v^i(N): ", round(m_cost, 2)) else ""
      .plot_mcstp(C_irred, arcs_N, results$value,
                  main_title = tit, sub_title = sub_tit)
    }
  }

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$value, names(results$value))
  vals <- results$value[ord]; noms <- names(results$value)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    coalitions = coals_names,
    C_irred = C_irred,
    v_i = values,
    contributions = results$contributions,
    shapley = results$value,
    total = m_cost,
    nperms = nrow(results$contributions),
    percentage = round((results$value / m_cost) * 100, 2),
    ranking = rank_star
  )

  class(output) <- "mcstp_irred"

  return(output)
}

#' @export
print.mcstp_irred <- function(x, ...) {
  v_S <- data.frame(
    S = x$coalitions,
    "v^i(S)" = x$v_i,
    row.names = NULL,
    check.names = FALSE
  )
  print(v_S)
  cat("\nSolution: Shapley value\n")
  print(round(x$shapley, 2))
  invisible(x)
}



#### Optimistic Game ####

#' Cost Allocations based on the Optimistic Game for MCSTP
#'
#' @description
#' This function computes the cost allocation associated with the optimistic game
#' (Bergantiños and Vidal-Puga, 2007b) for a minimum cost spanning tree problem
#' \eqn{(N_0, C)}. The optimistic game, denoted as \eqn{(N, v_C^o)} \eqn{(}or simply
#' \eqn{(N, v^o)}\eqn{)}, evaluates the cost of each coalition \eqn{S \subseteq N}
#' assuming that agents in \eqn{N \setminus S} are already connected. Thus,
#' agents in \eqn{S} can connect to the source through agents in \eqn{N \setminus S}
#' for free. When the Shapley value is applied as the solution concept to this game,
#' the resulting allocation coincides with the Shapley value of the irreducible game \eqn{(N, v_C^i)}
#' (Bergantiños and Vidal-Puga, 2007b).

#'
#' @details
#' Let \eqn{T := N \setminus S}. The optimistic game \eqn{(N, v_C^o)} is a
#' cooperative game with transferable utility (TU game) defined for each
#' coalition \eqn{S \subseteq N} as:
#'
#' \deqn{v_C^o(S) = m(S_0, C^T),}
#'
#' where \eqn{m(S_0, C^T)} is the cost of the minimum spanning tree of the
#' subproblem induced by the agents in \eqn{S} and the source, using the
#' modified cost matrix \eqn{C^T}. In this matrix, the cost between agents
#' in \eqn{S} remains the same, i.e., \eqn{c_{ij}^T = c_{ij}} for all \eqn{i,j \in S},
#' but the cost to connect any agent \eqn{i \in S} to the source is updated as:
#'
#' \deqn{c_{0i}^T = \min_{j \in T \cup \{0\}} c_{ji}.}
#'
#' Applying the Shapley value to this TU game assigns to each agent \eqn{i \in N}
#' the average of their marginal contributions over all possible permutations
#' \eqn{\pi \in \Pi_N}:
#'
#' \deqn{Sh_i(N, v^o) = \displaystyle{\frac{1}{n!} \sum_{\pi \in \Pi_N} [v^o(Pre(i, \pi) \cup \{i\}) - v^o(Pre(i, \pi))]},}
#'
#' where \eqn{Pre(i, \pi)} is the set of players that precede agent \eqn{i}
#' in the ordering \eqn{\pi}.
#'
#' @inheritParams private_game
#'
#' @note
#' The optimistic game associated with any MCST problem coincides with the optimistic
#' game associated with its irreducible form, i.e.,
#'
#' \deqn{v_C^o = v_{\bar{C}}^o.}
#'
#' If \eqn{(N_0, \bar{C})} is irreducible, the private and optimistic games, \eqn{v^p}
#' and \eqn{v^o}, are dual: \eqn{v^p(S) + v^o(N \setminus S) = m(N_0, C)} for all \eqn{S \subset N}.
#'
#' Furthermore, as \eqn{Sh(N, v_C^o) = Sh(N, v_C^i)}, applying the
#' Shapley value to \eqn{v^o} is simply another way of obtaining the
#' \code{\link{folk_rule}} (Bergantiños and Vidal-Puga, 2007b).
#'
#' @return A list containing:
#' \itemize{
#'    \item \code{coalitions}: a character vector representing all possible coalitions \eqn{S \subseteq N} ordered by cardinality.
#'    \item \code{C_opt}: a list of the optimistic cost matrices \eqn{C^T} for each coalition \eqn{S \subseteq N}.
#'    \item \code{v_o}: a vector containing the characteristic function values \eqn{v^o(S)} for all coalitions.
#'    \item \code{contributions}: a data frame with the marginal contributions for all permutations \eqn{\pi \in \Pi_N}.
#'    \item \code{shapley}: the Shapley value allocation vector.
#'    \item \code{total}: the total cost of the MCST, \eqn{v^o(N)}.
#'    \item \code{nperms}: the number of permutations \eqn{n!} used to compute marginal contributions and the resulting Shapley value.
#'    \item \code{percentage}: the share of the total cost allocated to each agent.
#'    \item \code{ranking}: ranking of agents by cost (from highest to lowest; ties marked with *).
#' }
#'
#' @seealso
#' \code{\link{private_game}} for a game based on a pessimistic approach.
#'
#' \code{\link{irred_game}} for an equivalent game under the Shapley value.
#'
#' \code{\link{alloc_rules}} for an overview of the available rules in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2007b) The optimistic TU game in minimum cost
#' spanning tree problems. Int J Game Theory 36(2):223–239
#'
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' @examples
#' # Simple vector input
#' opt_game(c(12, 15, 20, 4, 6, 8), sol = "shapley", draw = TRUE)
#'
#' # Input with infinite costs (disconnected nodes)
#' C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
#'            Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
#' opt_game(C_inf, sol = "shapley")
#'
#' # Matrix input
#' C_mat <- matrix(c(0, 12, 15, 12,
#'                  12,  0,  4,  6,
#'                  15,  4,  0,  8,
#'                  12,  6,  8,  0), byrow = TRUE, ncol = 4)
#' opt_game(C_mat, draw = TRUE, which = 2, titles = FALSE)
#'
#' @concept Cooperative Games
#' @concept MCSTP
#'
#' @export

opt_game <- function(C, sol = "shapley", draw = FALSE, which = c(1, 2), titles = TRUE) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1  # Number of agents n = |N|

  # Build the characteristic function v^o(S) using binary indexing
  ncoals <- 2^n - 1
  vvals <- numeric(ncoals)

  C_T <- list() # List to store C^{N \ S} for each coalition S

  for (i in 1:ncoals) {
    # Decode binary index i to find members of coalition S
    S <- which(as.logical(intToBits(i)[1:n]))

    # Identify the complementary coalition T = N \ S
    T_agents <- setdiff(1:n, S)

    sub_idx <- c(1, S + 1) # Include the source (0)
    T0_idx <- c(1, T_agents + 1) # T_0

    # Standard matrix for S
    sub_C <- C_mat[sub_idx, sub_idx, drop = FALSE]

    # Modify the connection to the source for each agent in S
    # c^T_{0i} = min_{j \in T_0} c_{ji}
    for (k in 1:length(S)) {
      s_idx <- S[k] + 1
      min_cost <- min(C_mat[s_idx, T0_idx])

      # Update the sub_C matrix (symmetric)
      sub_C[1, k + 1] <- min_cost
      sub_C[k + 1, 1] <- min_cost
    }

    # Store the modified matrix C^T for current coalition
    C_T[[i]] <- sub_C

    # Compute MST cost for this specific subset
    vvals[i] <- .get_cost(sub_C)
  }


  ## Execution ##

  # Calculate the allocation based on the chosen solution concept
  if (sol == "shapley") {
    results <- .compute_shapley(vvals, n)
  } else {
    stop("Solution not implemented. Currently, only 'shapley' is available")
  }

  m_cost <- vvals[ncoals] # Total cost v^o(N)

  # Generate all possible coalitions ordered by cardinality for final output
  coalitions <- unlist(lapply(1:n, function(k) combn(n, k, simplify = FALSE)), recursive = FALSE)
  coals_names <- sapply(coalitions, function(S) paste0("{", paste(S, collapse = ","), "}"))
  coals_names[length(coals_names)] <- "N"
  idx <- sapply(coalitions, function(S) sum(2^(S - 1)))

  # To display characteristic function v^o(S)
  values <- vvals[idx]
  # v_S <- as.data.frame(t(values))
  # colnames(v_S) <- coals_names; rownames(v_S) <- "v^o(S)"
  # v_S <- data.frame(
  #   S = coals_names,
  #   "v^o(S)" = values,
  #   row.names = NULL,
  #   check.names = FALSE
  # )

  # To display matrices C^T
  C_opt <- list()
  for (S in coalitions) {
    S_names <- paste0("S_", paste(S, collapse = ""))
    C_opt[[S_names]] <- C_T[[sum(2^(S - 1))]]
  }

  # cat("Optimistic Game\n")
  # cat("Characteristic function v^o(S):\n")
  # print(v_S)
  # cat("---------------\n")
  # cat(paste0("Total Cost v^o(N): ", round(m_cost, 2), "\n"))
  # cat("\nSolution: Shapley value\n")
  # print(round(results$value, 2))

  # if (n <= 3) {
  #   cat("\nDisplaying marginal contributions:\n")
  #   nperms <- nrow(results$contributions)
  #   cat(paste0("Computed ", nperms, " permutations\n"))
  #   print(results$contributions)
  # } else {
  #   nperms <- factorial(n)
  #   cat(paste0("\nMarginal contributions: computed ", nperms, " permutations (too many to display individually)\n"))
  # }


  ## Visualization ##

  if (draw) {

    if (is.null(which)) which <- c(1, 2)
    show_details <- 1 %in% which
    show_main <- 2 %in% which

    if (n > 1 && n <= 3 && show_details) {

      rows <- if(n == 3) 2 else 1
      cols <- if(n == 3) 3 else 2
      old_par <- par(mfrow = c(rows, cols), oma = c(0, 0, if(titles) 5 else 0, 0))
      on.exit(par(old_par), add = TRUE)

      for (i in 1:(ncoals - 1)) {
        S <- coalitions[[i]]

        # Recalculate sub_C for visualization with modified source connections
        T_agents <- setdiff(1:n, S)
        sub_idx <- c(1, S + 1)
        T0_idx <- c(1, T_agents + 1)

        C_inf <- matrix(Inf, nrow = n + 1, ncol = n + 1, dimnames = list(0:n, 0:n))

        sub_C <- C_mat[sub_idx, sub_idx, drop = FALSE]

        for (k in 1:length(S)) {
          s_idx <- S[k] + 1
          min_cost <- min(C_mat[s_idx, T0_idx])
          sub_C[1, k + 1] <- min_cost
          sub_C[k + 1, 1] <- min_cost
        }

        C_inf[sub_idx, sub_idx] <- sub_C
        arcs_S <- .get_arcs(sub_C)

        .plot_mcstp(C_inf, arcs_S, NULL,
                    main_title = "", sub_title = paste0("S = {", paste(S, collapse=","), "}")) #   |  v^o(S): ", round(values[i], 2)
      }

      if (titles) {
        mtext("MCST for each subcoalition S",
              side = 3, line = 1.5, cex = 1.2, font = 2, outer = TRUE)
        mtext("Optimistic Game (N, v^o)", side = 3, line = 0.1, cex = 0.7, family = "sans", font = 3,
              col = "#444444", outer = TRUE) # paste0("Irreducible Game  |  Total Network Cost v^i(N): ", round(m_cost, 2))
      }

      par(old_par)
    }

    if (n > 1 && n <= 3 && show_details && show_main) {
      devAskNewPage(TRUE)
      on.exit(devAskNewPage(FALSE), add = TRUE)
    }

    if (show_main) {
      arcs_N <- .get_arcs(C_mat)
      tit <- if(titles) "Optimistic Game (N, v^o)" else ""
      sub_tit <- if(titles) paste0("Solution: Shapley value  |  Total Network Cost v^o(N): ", round(m_cost, 2)) else ""
      .plot_mcstp(C_mat, arcs_N, results$value,
                  main_title = tit, sub_title = sub_tit)
    }
  }

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$value, names(results$value))
  vals <- results$value[ord]; noms <- names(results$value)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    coalitions = coals_names,
    C_opt = C_opt,
    v_o = values,
    contributions = results$contributions,
    shapley = results$value,
    total = m_cost,
    nperms = nrow(results$contributions),
    percentage = round((results$value / m_cost) * 100, 2),
    ranking = rank_star
  )

  class(output) <- "mcstp_opt"

  return(output)
}

#' @export
print.mcstp_opt <- function(x, ...) {
  v_S <- data.frame(
    S = x$coalitions,
    "v^o(S)" = x$v_o,
    row.names = NULL,
    check.names = FALSE
  )
  print(v_S)
  cat("\nSolution: Shapley value\n")
  print(round(x$shapley, 2))
  invisible(x)
}



#### Public Game ####

#' Cost Allocations based on the Public Game for MCSTP
#'
#' @description
#' This function computes the cost allocation associated with the public game
#' (Bogomolnaia and Moulin, 2010) for a minimum cost spanning tree problem
#' \eqn{(N_0, C)}. The public game, denoted as \eqn{(N, v_C^u)} \eqn{(}or simply
#' \eqn{(N, v^u)}\eqn{)}, evaluates the cost of each coalition \eqn{S \subseteq N}
#' assuming that there are no property rights on nodes. Therefore, agents in
#' \eqn{S} can use the nodes of agents in \eqn{N \setminus S} to connect to
#' the source, paying the corresponding arc costs.
#'
#' @details
#' The public game \eqn{(N, v_C^u)} is a cooperative game with transferable
#' utility (TU game) where the cost for each coalition \eqn{S \subseteq N}
#' is defined as the minimum private game cost among all its supersets
#' \eqn{T \supseteq S}, i.e.:
#'
#' \deqn{v_C^u(S) = \min_{S \subseteq T} v_C^p(T),}
#'
#' where \eqn{v_C^p(\cdot)} is the characteristic function of the private game.
#' This means computing \eqn{v^u(S)} requires calculating the private game costs
#' for all possible supersets of \eqn{S}.
#'
#' Applying the Shapley value to this TU game assigns to each agent \eqn{i \in N}
#' the average of their marginal contributions over all possible permutations
#' \eqn{\pi \in \Pi_N}:
#'
#' \deqn{Sh_i(N, v^u) = \displaystyle{\frac{1}{n!} \sum_{\pi \in \Pi_N} [v^u(Pre(i, \pi) \cup \{i\}) - v^u(Pre(i, \pi))]},}
#'
#' where \eqn{Pre(i, \pi)} is the set of players that precede agent \eqn{i}
#' in the ordering \eqn{\pi}.
#'
#' @inheritParams private_game
#'
#' @note
#' The core of \eqn{v^u} coincides with the core of \eqn{v^p} restricted
#' to nonnegative cost shares. Consequently, classical rules such as Bird,
#' Dutta-Kar, or the folk rule belong to the core of the public game
#' for any MCST problem.
#'
#' However, unlike in the private game, the Shapley value of the public game
#' does not necessarily belong to the core of \eqn{v^u} (Trudeau and Vidal-Puga, 2019).
#'
#' @return A list containing:
#' \itemize{
#'    \item \code{coalitions}: a character vector representing all possible coalitions \eqn{S \subseteq N} ordered by cardinality.
#'    \item \code{v_u}: a vector containing the characteristic function values \eqn{v^u(S)} for all coalitions.
#'    \item \code{v_p}: a vector containing the characteristic function values \eqn{v^p(S)} for all coalitions.
#'    \item \code{best_T}: a data frame mapping each coalition \eqn{S} to the optimal superset \eqn{T} that minimizes its cost.
#'    \item \code{contributions}: a data frame with the marginal contributions for all permutations \eqn{\pi \in \Pi_N}.
#'    \item \code{shapley}: the Shapley value allocation vector.
#'    \item \code{total}: the total cost of the network, \eqn{v^u(N)}.
#'    \item \code{nperms}: the number of permutations \eqn{n!} used to compute marginal contributions and the resulting Shapley value.
#'    \item \code{percentage}: the share of the total cost allocated to each agent.
#'    \item \code{ranking}: ranking of agents by cost (from highest to lowest; ties marked with *).
#' }
#'
#' @seealso
#' \code{\link{bird_rule}}, \code{\link{dk_rule}}, \code{\link{folk_rule}}
#' for classical rules belonging to the core of \eqn{v^u}.
#'
#' \code{\link{private_game}} for the pessimistic alternative to the public game.
#'
#' \code{\link{irred_game}} for a game that can be defined
#' from either the private or the public approach.
#'
#' \code{\link{alloc_rules}} for an overview of the available rules in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' Bogomolnaia A, Moulin H (2010) Sharing a minimal cost spanning tree:
#' beyond the folk solution. Games and Econom Behav 69(2):238–248
#'
#' Trudeau C, Vidal-Puga J (2019) The Shapley value in minimum cost
#' spanning tree problems. In: Algaba E, Fragnelli V, Sánchez-Soriano
#' J (eds) Chapters in Game Theory: The Shapley value, chapter 24.
#' CRC Press, Taylor & Francis Group
#'
#' @examples
#' # Simple vector input
#' public_game(c(12, 15, 20, 4, 6, 8), sol = "shapley", draw = TRUE)
#'
#' # Input with infinite costs (disconnected nodes)
#' C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
#'            Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
#' public_game(C_inf, sol = "shapley")
#'
#' # Matrix input
#' C_mat <- matrix(c(0, 12, 15, 12,
#'                  12,  0,  4,  6,
#'                  15,  4,  0,  8,
#'                  12,  6,  8,  0), byrow = TRUE, ncol = 4)
#' public_game(C_mat, draw = TRUE, which = 2, titles = FALSE)
#'
#' @concept Cooperative Games
#' @concept MCSTP
#'
#' @export

public_game <- function(C, sol = "shapley", draw = FALSE, which = c(1, 2), titles = TRUE) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1  # Number of agents n = |N|

  # Compute the private game values v^p(T) for all coalitions
  ncoals <- 2^n - 1
  vvals_p <- numeric(ncoals)

  for (i in 1:ncoals) {
    # Decode binary index i to find members of T
    T_set <- which(as.logical(intToBits(i)[1:n]))
    sub_idx <- c(1, T_set + 1) # Include the source (0)
    sub_C <- C_mat[sub_idx, sub_idx, drop = FALSE]
    vvals_p[i] <- .get_cost(sub_C)
  }

  # Compute the public game values v^u(S)
  vvals <- numeric(ncoals)
  T_opt <- integer(ncoals) # To store the optimal superset T for drawing

  for (i in 1:ncoals) {
    # Find all supersets T of S
    supersets <- which(bitwAnd(1:ncoals, i) == i)

    # v^u(S) is the minimum v^p(T) among all supersets T
    min_idx <- supersets[which.min(vvals_p[supersets])]
    vvals[i] <- vvals_p[min_idx]
    T_opt[i] <- min_idx
  }


  ## Execution ##

  # Calculate the allocation based on the chosen solution concept
  if (sol == "shapley") {
    results <- .compute_shapley(vvals, n)
  } else {
    stop("Solution not implemented. Currently, only 'shapley' is available")
  }

  m_cost <- vvals[ncoals] # Total cost v^u(N)

  # Generate all possible coalitions ordered by cardinality for final output
  coalitions <- unlist(lapply(1:n, function(k) combn(n, k, simplify = FALSE)), recursive = FALSE)
  coals_names <- sapply(coalitions, function(S) paste0("{", paste(S, collapse = ","), "}"))
  coals_names[length(coals_names)] <- "N"
  idx <- sapply(coalitions, function(S) sum(2^(S - 1)))

  # Characteristic function v^p(S) for final output
  values_p <- vvals_p[idx]

  # Optimal superset for each coalition S for final output
  T_names <- character(2^n - 1)
  T_names[idx] <- coals_names

  # arg min_{ScT} v^p(T) for final output
  argmin_T <- T_names[T_opt[idx]]
  # names(argmin_T) <- coals_names
  # best_T <- noquote(argmin_T)
  best_T <- data.frame(
    S = coals_names,
    T = as.character(argmin_T),
    row.names = NULL,
    check.names = FALSE
  )

  # To display characteristic function v^u(S)
  values <- vvals[idx]
  # v_S <- as.data.frame(t(values))
  # colnames(v_S) <- coals_names; rownames(v_S) <- "v^u(S)"
  # v_S <- data.frame(
  #   S = coals_names,
  #   "v^u(S)" = values,
  #   row.names = NULL,
  #   check.names = FALSE
  # )

  # cat("Public Game\n")
  # cat("Characteristic function v^u(S):\n")
  # print(v_S)
  # cat("---------------\n")
  # cat(paste0("Total Cost v^u(N): ", round(m_cost, 2), "\n"))
  # cat("\nSolution: Shapley value\n")
  # print(round(results$value, 2))

  # if (n <= 3) {
  #   cat("\nDisplaying marginal contributions:\n")
  #   nperms <- nrow(results$contributions)
  #   cat(paste0("Computed ", nperms, " permutations\n"))
  #   print(results$contributions)
  # } else {
  #   nperms <- factorial(n)
  #   cat(paste0("\nMarginal contributions: computed ", nperms, " permutations (too many to display individually)\n"))
  # }


  ## Visualization ##

  if (draw) {

    if (is.null(which)) which <- c(1, 2)
    show_details <- 1 %in% which
    show_main <- 2 %in% which

    if (n > 1 && n <= 3 && show_details) {

      rows <- if(n == 3) 2 else 1
      cols <- if(n == 3) 3 else 2
      old_par <- par(mfrow = c(rows, cols), oma = c(0, 0, if(titles) 5 else 0, 0))
      on.exit(par(old_par), add = TRUE)

      for (i in 1:(ncoals - 1)) {
        S <- coalitions[[i]]

        # We retrieve the optimal superset T that provided the minimum cost
        T_set <- which(as.logical(intToBits(T_opt[sum(2^(S - 1))])[1:n]))
        sub_idx <- c(1, T_set + 1)

        C_inf <- matrix(Inf, nrow = n + 1, ncol = n + 1, dimnames = list(0:n, 0:n))

        C_inf[sub_idx, sub_idx] <- C_mat[sub_idx, sub_idx]
        sub_C <- C_mat[sub_idx, sub_idx, drop = FALSE]
        arcs_T <- .get_arcs(sub_C)

        .plot_mcstp(C_inf, arcs_T, NULL,
                    main_title = "", sub_title = paste0("S = {", paste(S, collapse=","), "}")) #   |  v^u(S): ", round(values[i], 2)
      }

      if (titles) {
        mtext("MCST for each subcoalition S",
              side = 3, line = 1.5, cex = 1.2, font = 2, outer = TRUE)
        mtext("Public Game (N, v^u)", side = 3, line = 0.1, cex = 0.7, family = "sans", font = 3,
              col = "#444444", outer = TRUE) # paste0("Public Game  |  Total Network Cost v^u(N): ", round(m_cost, 2))
      }

      par(old_par)
    }

    if (n > 1 && n <= 3 && show_details && show_main) {
      devAskNewPage(TRUE)
      on.exit(devAskNewPage(FALSE), add = TRUE)
    }

    if (show_main) {
      arcs_N <- .get_arcs(C_mat)
      tit <- if(titles) "Public Game (N, v^u)" else ""
      sub_tit <- if(titles) paste0("Solution: Shapley value  |  Total Network Cost v^u(N): ", round(m_cost, 2)) else ""
      .plot_mcstp(C_mat, arcs_N, results$value,
                  main_title = tit, sub_title = sub_tit)
    }
  }

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$value, names(results$value))
  vals <- results$value[ord]; noms <- names(results$value)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    coalitions = coals_names,
    v_u = values,
    v_p = values_p,
    best_T = best_T,
    contributions = results$contributions,
    shapley = results$value,
    total = m_cost,
    nperms = nrow(results$contributions),
    percentage = round((results$value / m_cost) * 100, 2),
    ranking = rank_star
  )

  class(output) <- "mcstp_public"

  return(output)
}

#' @export
print.mcstp_public <- function(x, ...) {
  v_S <- data.frame(
    S = x$coalitions,
    "v^u(S)" = x$v_u,
    row.names = NULL,
    check.names = FALSE
  )
  print(v_S)
  cat("\nSolution: Shapley value\n")
  print(round(x$shapley, 2))
  invisible(x)
}



#### Cycle-complete Game ####

#' Cost Allocations based on the Cycle-complete Game for MCSTP
#'
#' @description
#' This function computes the cost allocation associated with the cycle-complete game
#' (Trudeau, 2012) for a minimum cost spanning tree problem \eqn{(N_0, C)}. The
#' cycle-complete game, denoted as \eqn{(N, v_C^c)} \eqn{(}or simply \eqn{(N, v^c)}\eqn{)}, is
#' defined as the private game associated with the cycle-complete network
#' \eqn{(N_0, C^*)}. When the Shapley value is applied as the solution
#' concept to this game, the resulting allocation is called \emph{cycle-complete solution}.
#'
#' @details
#' The cycle-complete game \eqn{(N, v_C^c)} is a cooperative game with
#' transferable utility (TU game) defined for each coalition \eqn{S \subseteq N} as:
#'
#' \deqn{v_C^c(S) = m(S_0, C^*),}
#'
#' where \eqn{m(S_0, C^*)} is the cost of the minimum spanning tree of the
#' subproblem induced by the agents in \eqn{S} and the source, evaluated over
#' the cycle-complete costs \eqn{C^*}. The cycle-complete cost is defined as
#' the minimum between the direct cost \eqn{c_{ij}} and the smallest value of
#' the maximum arc cost over all cycles containing \eqn{i} and \eqn{j}, i.e.,
#' \eqn{c_{ij}^* = \min(c_{ij}, \min_{f \in \mathcal{C}(i,j)} \max_{e \in f} c_e),}
#' where \eqn{\mathcal{C}(i,j)} denotes the set of all cycles that include both \eqn{i} and \eqn{j}.
#'
#' For computational efficiency, \eqn{C^*} is obtained via
#' the irreducible matrices \eqn{\bar{C}} of the subgraphs as:
#'
#' \deqn{c_{ij}^* = \max_{k \in N \setminus \{i,j\}} {\bar{c}}_{ij}^{N \setminus \{k\}} \text{ for } i,j \in N,}
#'
#' and
#'
#' \deqn{c_{0i}^* = \max_{k \in N \setminus \{i\}} {\bar{c}}_{0i}^{N \setminus \{k\}} \text{ for } i \in N,}
#'
#' where \eqn{{\bar{C}}^{N \setminus \{k\}}} is the irreducible cost matrix of
#' the network excluding agent \eqn{k}. If no agents can be removed, the
#' irreducible cost of the full network is used.
#'
#' Applying the Shapley value to this TU game, known as the
#' \emph{cycle-complete solution}, assigns to each agent \eqn{i \in N}
#' the average of their marginal contributions over all possible permutations
#' \eqn{\pi \in \Pi_N}:
#'
#' \deqn{Sh_i(N, v^c) = \displaystyle{\frac{1}{n!} \sum_{\pi \in \Pi_N} [v^c(Pre(i, \pi) \cup \{i\}) - v^c(Pre(i, \pi))]},}
#'
#' where \eqn{Pre(i, \pi)} is the set of players that precede agent \eqn{i}
#' in the ordering \eqn{\pi}.
#'
#' @inheritParams private_game
#'
#' @note
#' The cycle-complete game is concave, which implies that its Shapley value
#' (the \emph{cycle-complete solution}) belongs to its core. Moreover,
#' \eqn{\mathrm{core}(N,v^c) \subseteq \mathrm{core}(N,v^p)}.
#'
#' @return A list containing:
#' \itemize{
#'    \item \code{coalitions}: a character vector representing all possible coalitions \eqn{S \subseteq N} ordered by cardinality.
#'    \item \code{C_cc}: the cycle-complete cost matrix \eqn{C^*}.
#'    \item \code{v_c}: a vector containing the characteristic function values \eqn{v^c(S)} for all coalitions.
#'    \item \code{contributions}: a data frame with the marginal contributions for all permutations \eqn{\pi \in \Pi_N}.
#'    \item \code{shapley}: the Shapley value allocation vector (\emph{cycle-complete solution}).
#'    \item \code{total}: the total cost of the network, \eqn{v^c(N)}.
#'    \item \code{nperms}: the number of permutations \eqn{n!} used to compute marginal contributions and the resulting Shapley value.
#'    \item \code{percentage}: the share of the total cost allocated to each agent.
#'    \item \code{ranking}: ranking of agents by cost (from highest to lowest; ties marked with *).
#' }
#'
#' @seealso
#' \code{\link{private_game}} for the pessimistic game based on the original costs.
#'
#' \code{\link{irred_game}} for the pessimistic game based on the irreducible cost matrix \eqn{\bar{C}}.
#'
#' \code{\link{alloc_rules}} for an overview of the available rules in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' Trudeau C (2012) A new stable and more responsible cost sharing solution
#' for mcst problems. Games Econom Behav 75(1):402–412
#'
#' @examples
#' # Simple vector input
#' cc_game(c(12, 15, 20, 4, 6, 8), sol = "shapley", draw = TRUE)
#'
#' # Input with infinite costs (disconnected nodes)
#' C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
#'            Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
#' cc_game(C_inf, sol = "shapley")
#'
#' # Matrix input
#' C_mat <- matrix(c(0, 12, 15, 12,
#'                  12,  0,  4,  6,
#'                  15,  4,  0,  8,
#'                  12,  6,  8,  0), byrow = TRUE, ncol = 4)
#' cc_game(C_mat, draw = TRUE, which = 2, titles = FALSE)
#'
#' @concept Cooperative Games
#' @concept MCSTP
#'
#' @export

cc_game <- function(C, sol = "shapley", draw = FALSE, which = c(1, 2), titles = TRUE) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  n <- nrow(C_mat) - 1  # Number of agents n = |N|

  # Compute the cycle-complete matrix C*
  C_cc <- .get_cycle_complete(C_mat)

  # Build the characteristic function v^c(S) using binary indexing over C*
  ncoals <- 2^n - 1
  vvals <- numeric(ncoals)

  for (i in 1:ncoals) {
    # Decode binary index i to find members of coalition S
    S <- which(as.logical(intToBits(i)[1:n]))
    sub_idx <- c(1, S + 1) # Include the source (node 0)
    sub_C <- C_cc[sub_idx, sub_idx, drop = FALSE]

    # Compute MST cost for this specific subset
    vvals[i] <- .get_cost(sub_C)
  }


  ## Execution ##

  # Calculate the allocation based on the chosen solution concept
  if (sol == "shapley") {
    results <- .compute_shapley(vvals, n)
  } else {
    stop("Solution not implemented. Currently, only 'shapley' is available")
  }

  m_cost <- vvals[ncoals] # Total cost v^c(N)

  # Generate all possible coalitions ordered by cardinality for final output
  coalitions <- unlist(lapply(1:n, function(k) combn(n, k, simplify = FALSE)), recursive = FALSE)
  coals_names <- sapply(coalitions, function(S) paste0("{", paste(S, collapse = ","), "}"))
  coals_names[length(coals_names)] <- "N"
  idx <- sapply(coalitions, function(S) sum(2^(S - 1)))

  # To display characteristic function v^c(S)
  values <- vvals[idx]
  # v_S <- as.data.frame(t(values))
  # colnames(v_S) <- coals_names; rownames(v_S) <- "v^c(S)"
  # v_S <- data.frame(
  #   S = coals_names,
  #   "v^c(S)" = values,
  #   row.names = NULL,
  #   check.names = FALSE
  # )

  # cat("Cycle-complete Game\n")
  # cat("Characteristic function v^c(S):\n")
  # print(v_S)
  # cat("---------------\n")
  # cat(paste0("Total Cost v^c(N): ", round(m_cost, 2), "\n"))
  # cat("\nSolution: Shapley value (cycle-complete solution)\n")
  # print(round(results$value, 2))

  # if (n <= 3) {
  #   cat("\nDisplaying marginal contributions:\n")
  #   nperms <- nrow(results$contributions)
  #   cat(paste0("Computed ", nperms, " permutations\n"))
  #   print(results$contributions)
  # } else {
  #   nperms <- factorial(n)
  #   cat(paste0("\nMarginal contributions: computed ", nperms, " permutations (too many to display individually)\n"))
  # }


  ## Visualization ##

  if (draw) {

    if (is.null(which)) which <- c(1, 2)
    show_details <- 1 %in% which
    show_main <- 2 %in% which

    if (n > 1 && n <= 3 && show_details) {

      rows <- if(n == 3) 2 else 1
      cols <- if(n == 3) 3 else 2
      old_par <- par(mfrow = c(rows, cols), oma = c(0, 0, if(titles) 5 else 0, 0))
      on.exit(par(old_par), add = TRUE)

      for (i in 1:(ncoals - 1)) {
        S <- coalitions[[i]]
        sub_idx <- c(1, S + 1)

        C_inf <- matrix(Inf, nrow = n + 1, ncol = n + 1, dimnames = list(0:n, 0:n))

        C_inf[sub_idx, sub_idx] <- C_cc[sub_idx, sub_idx]
        sub_C <- C_cc[sub_idx, sub_idx, drop = FALSE]
        arcs_S <- .get_arcs(sub_C)

        .plot_mcstp(C_inf, arcs_S, NULL,
                    main_title = "", sub_title = paste0("S = {", paste(S, collapse=","), "}")) #   |  v^c(S): ", round(values[i], 2)
      }

      if (titles) {
        mtext("MCST for each subcoalition S",
              side = 3, line = 1.5, cex = 1.2, font = 2, outer = TRUE)
        mtext("Cycle-complete Game (N, v^c)", side = 3, line = 0.1, cex = 0.7, family = "sans", font = 3,
              col = "#444444", outer = TRUE) # paste0("Cycle-complete Game  |  Total Network Cost v^c(N): ", round(m_cost, 2))
      }

      par(old_par)
    }

    if (n > 1 && n <= 3 && show_details && show_main) {
      devAskNewPage(TRUE)
      on.exit(devAskNewPage(FALSE), add = TRUE)
    }

    if (show_main) {
      arcs_N <- .get_arcs(C_cc)
      tit <- if(titles) "Cycle-complete Game (N, v^c)" else ""
      sub_tit <- if(titles) paste0("Cycle-complete solution  |  Total Network Cost v^c(N): ", round(m_cost, 2)) else ""
      .plot_mcstp(C_cc, arcs_N, results$value,
                  main_title = tit, sub_title = sub_tit)
    }
  }

  # Generates a ranking based on allocations, using stars (*) for ties
  ord <- order(-results$value, names(results$value))
  vals <- results$value[ord]; noms <- names(results$value)[ord]
  tab <- table(vals)
  rep_vals <- sort(as.numeric(names(tab[tab > 1])), decreasing = TRUE)

  rank_star <- noquote(sapply(1:length(vals), function(i) {
    if (vals[i] %in% rep_vals) {
      paste0(noms[i], strrep("*", which(rep_vals == vals[i])))
    } else noms[i]
  }))


  ## Output ##

  output <- list(
    coalitions = coals_names,
    C_cc = C_cc,
    v_c = values,
    contributions = results$contributions,
    shapley = results$value,
    total = m_cost,
    nperms = nrow(results$contributions),
    percentage = round((results$value / m_cost) * 100, 2),
    ranking = rank_star
  )

  class(output) <- "mcstp_cc"

  return(output)
}

#' @export
print.mcstp_cc <- function(x, ...) {
  v_S <- data.frame(
    S = x$coalitions,
    "v^c(S)" = x$v_c,
    row.names = NULL,
    check.names = FALSE
  )
  print(v_S)
  cat("\nSolution: Shapley value\n")
  print(round(x$shapley, 2))
  invisible(x)
}



#### Internal functions ####

# Prepares the cost matrix from a vector or matrix, ensuring row/column names
.prepare_matrix <- function(x) {
  if (is.vector(x)) {
    # Solve for N in the equation length(x) = N * (N - 1) / 2
    N <- (1 + sqrt(1 + 8 * length(x))) / 2
    m <- matrix(0, N, N)
    m[lower.tri(m)] <- x
    m <- m + t(m)
    rownames(m) <- colnames(m) <- as.character(0:(N-1))
    return(m)
  }
  if (is.matrix(x)) {
    if (is.null(rownames(x))) {
      rownames(x) <- colnames(x) <- as.character(0:(nrow(x)-1))
    }
    return(x)
  }
  stop("Input must be a numeric vector or a square matrix.")
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


# Calculation of the Shapley value via marginal contributions over permutations
.compute_shapley <- function(vvals, n) {

  # Generate all possible permutations of the n players
  perms <- gtools::permutations(n, n, 1:n)
  nperms <- nrow(perms)

  # Initialize matrix for marginal contributions
  contributions <- matrix(0, nrow = nperms, ncol = n, dimnames = list(NULL, as.character(1:n)))

  # Inner helper to retrieve v(S) from the binary-indexed vector
  .get_v <- function(S) {
    if (length(S) == 0) return(0)
    # Map subset S to its binary position
    idx <- sum(2^(S - 1))
    return(vvals[idx])
  }

  # Loop through each permutation to calculate individual contributions
  for (p in 1:nperms) {
    pi <- perms[p, ]
    predecessors <- c()

    for (i in pi) {
      # Marginal contribution: v(Pre(i,pi) U {i}) - v(Pre(i,pi))
      v_with_i <- .get_v(c(predecessors, i))
      v_without_i <- .get_v(predecessors)

      # Store the difference for player i in this specific permutation
      contributions[p, i] <- v_with_i - v_without_i

      # Update the set of players who have already arrived
      predecessors <- c(predecessors, i)
    }
  }

  # Shapley value: average of all marginal contributions
  sh <- colMeans(contributions)

  # Convert contributions matrix to data frame
  contrib_df <- data.frame(
    pi = apply(perms, 1, paste, collapse = "-"),
    contributions,
    row.names = NULL,
    check.names = FALSE
  )

  return(list(
    value = sh,
    contributions = contrib_df
  ))
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
.plot_mcstp <- function(adj, mst_arcs, allocation, main_title, sub_title) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' needed for drawing. Please install it.")
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
    # stage <- mst_arcs$stage[k]

    for(l in 1:nrow(arcs)) {
      if((arcs[l,1] == i && arcs[l,2] == j) ||
         (arcs[l,1] == j && arcs[l,2] == i)) {
        arc_colors[l] <- "red"
        arc_widths[l] <- 2
        # arc_labels[l] <- paste0(arc_labels[l], " [", stage, "]")
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


# Modified 29/04 13:23
