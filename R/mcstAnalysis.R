#### Analysis functions ####

  # Core Check:             mcstCoreCheck + plot
  # Core Points:            mcstCorePoints + plot
  # Core Plot:              mcstCorePlot
  # Stability Range:        mcstStabilityRange
  # Sensitivity:            mcstSensitivity + plot
  # Comparison              mcstCompare + plot


#### Core Check ####

#' Core Stability Check for Cost Allocations
#'
#' @description
#' This function evaluates whether a given cost allocation belongs to the \emph{core}
#' of a cooperative game. It checks for both efficiency and coalitional rationality,
#' identifying any blocking coalitions if the allocation is unstable.
#'
#' @details
#' An allocation \eqn{x \in \mathbb{R}^N} belongs to the \emph{core} of a cooperative game \eqn{(N, v)} if
#' it satisfies two fundamental conditions:
#'
#' \enumerate{
#' \item Efficiency: the total cost is completely allocated among the agents, i.e.,
#' \deqn{\sum_{i \in N} x_i = v(N).}
#'
#' \item Coalitional rationality: no coalition pays more than its stand-alone cost, i.e.,
#' \deqn{\sum_{i \in S} x_i \le v(S) \quad \forall S \subset N.}
#'
#' }
#'
#' If the coalitional rationality condition is violated for any coalition \eqn{S},
#' \eqn{S} is considered a \emph{blocking coalition}, as its members would have an incentive
#' to leave the grand coalition and form their own sub-network. The excess
#' \eqn{x(S) - v(S)} quantifies the magnitude of this violation.
#'
#' The function can be executed in two ways: by passing a complete \code{game} object
#' (which must have been computed with \code{method = "exact"}), or by manually providing
#' the proposed \code{allocation} vector and the characteristic function \code{v}.
#'
#' The \code{plot} method visualizes the excess of the blocking coalitions (if any)
#' using a bar chart.
#'
#' @param game list; an object containing a cooperative game computed with
#' \code{method = "exact"} (e.g., from \code{\link{mcstGamePrivate}}). Default is \code{NULL}.
#' @param allocation numeric vector; a proposed cost allocation for each agent.
#' Required if \code{game} is \code{NULL}.
#' @param v numeric vector; the characteristic function values \eqn{v(S)} for all
#' \eqn{2^n - 1} coalitions, ordered by cardinality. Required if \code{game} is \code{NULL}.
#' @param tol numeric; tolerance for floating-point comparisons to avoid precision
#' issues. Default is \code{1e-7}.
#'
#' @return A list containing:
#' \itemize{
#'    \item \code{in_core}: logical; if \code{TRUE}, the allocation is stable (belongs to the \emph{core}).
#'    \item \code{is_efficient}: logical; if \code{TRUE}, the sum of the allocation equals the total cost.
#'    \item \code{is_rational}: logical; if \code{TRUE}, no blocking coalitions exist.
#'    \item \code{blocking_coals}: a data frame detailing the blocking coalitions \eqn{S}, their allocated cost \eqn{x(S)}, their value \eqn{v(S)}, and the excess. Empty if the allocation is rational.
#'    \item \code{n}: the number of agents.
#'    \item \code{sum_x}: the sum of the allocated costs.
#'    \item \code{total}: the total cost of the game \eqn{v(N)}.
#' }
#'
#' @seealso
#' \code{\link{mcstGamePrivate}}, \code{\link{mcstGameIrred}} for computing cooperative games
#' that can be passed to this function.
#'
#' \code{\link{mcstRules}} for an overview of the available rules and analysis tools in the package.
#'
#' @references
#' Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100
#'
#' @examples
#' # Example 1: Using manual vectors
#' vvals <- c(10, 12, 15, 20, 22, 25, 30) # 3 agents: v(1), v(2), v(3), v(12), ...
#' # Stable check
#' alloc1 <- c(8, 10, 12)
#' check1 <- mcstCoreCheck(allocation = alloc1, v = vvals); check1
#' # Unstable check
#' alloc2 <- c(15, 10, 5)
#' check2 <- mcstCoreCheck(allocation = alloc2, v = vvals); check2
#' plot(check2)
#' \donttest{
#' # Example 2: Using a game object
#' C <- matrix(c(0, 12, 15, 12, 0, 4, 15, 4, 0), byrow = TRUE, ncol = 3)
#' # Check if the Shapley value of the private game is in the core
#' mcstCoreCheck(game = mcstGamePrivate(C))
#' }
#'
#' @concept Core Stability
#' @concept Cooperative Games
#'
#' @export

mcstCoreCheck <- function(game = NULL, allocation = NULL, v = NULL, tol = 1e-7) {

  # Check inputs
  if (!is.null(game) && is.list(game)) {
    if (is.null(game$method) || game$method != "exact") {
      stop("Required a 'game' object computed with method = 'exact'")
    }
    alloc <- game$allocation
    v_game <- names(game)[grep("^v_", names(game))][1]
    vvals <- game[[v_game]]
    m_cost <- game$total

  } else if (!is.null(allocation) && !is.null(v)) {
    alloc <- allocation
    vvals <- v
    m_cost <- vvals[length(vvals)]

  } else {
    stop("Provide either a 'game' object, or both 'allocation' and 'v' parameters")
  }


  ## Execution ##

  n <- length(alloc)

  # Generate coalitions in cardinality order
  coals <- unlist(lapply(1:n, function(k) combn(n, k, simplify = FALSE)), recursive = FALSE)
  ncoals <- 2^n - 1

  # Check efficiency: sum of allocations equals total cost
  sum_alloc <- sum(alloc)
  is_efficient <- abs(sum_alloc - m_cost) < tol

  # Check coalitional rationality: x(S) <= v(S)
  blocking_list <- lapply(1:ncoals, function(k) {
    S <- coals[[k]]
    x_S <- sum(alloc[S])
    v_S <- vvals[k]

    if (x_S > v_S + tol && length(S) < n) {
      data.frame(
        S = paste0("{", paste(S, collapse = ","), "}"),
        "x(S)" = round(x_S, 2),
        "v(S)" = round(v_S, 2),
        excess = round(x_S - v_S, 2),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })

  blocking_coals <- do.call(rbind, blocking_list)
  if (is.null(blocking_coals)) {
    blocking_coals <- data.frame(S = character(), "x(S)" = numeric(),
                                 "v(S)" = numeric(), excess = numeric(),
                                 check.names = FALSE)
  }

  is_rational <- nrow(blocking_coals) == 0


  ## Output ##

  output <- list(
    in_core = (is_efficient && is_rational),
    is_efficient = is_efficient,
    is_rational = is_rational,
    blocking_coals = blocking_coals,
    n = n,
    sum_x = sum_alloc,
    total = m_cost
  )

  class(output) <- "mcstp_core"

  return(output)
}

#' @export
print.mcstp_core <- function(x, ...) {
  cat("-------------------------\n")
  cat(" Core Stability Analysis\n")
  cat("-------------------------\n")
  cat("Players (n) :", x$n, "\n")

  eff_status <- if(x$is_efficient) "Passed" else "Failed"
  cat("Efficiency  :", eff_status, sprintf("(Sum x: %.2f | v(N): %.2f)\n", x$sum_x, x$total))

  rat_status <- if(x$is_rational) "Passed" else "Failed"
  cat("Rationality :", rat_status, "\n\n")

  if (x$in_core) {
    cat("Result: the allocation IS IN THE CORE (Stable)\n")
  } else {
    cat("Result: the allocation IS NOT IN THE CORE (Unstable)\n")
    if (nrow(x$blocking_coals) > 0) {
      cat("\nBlocking Coalitions (x(S) > v(S)):\n")
      print(x$blocking_coals, row.names = FALSE)
    }
  }
  invisible(x)
}

#' @rdname mcstCoreCheck
#' @param x the output from \code{mcstCoreCheck}.
#' @param col character string specifying the color for the excess bars. Default is \code{"firebrick"}.
#' @param titles logical; if \code{TRUE} (default), adds a main title and a y-axis label to the plot.
#' @param ... additional graphical parameters passed to \code{\link[graphics]{barplot}}.
#' @export
plot.mcstp_core <- function(x, col = "firebrick", titles = TRUE, ...) {
  if (x$is_rational) {
    message("The allocation is rational: no 'blocking coalitions' to plot")
    return(invisible(NULL))
  }

  plot_data <- x$blocking_coals[order(x$blocking_coals$excess, decreasing = TRUE), ]
  barplot(plot_data$excess,
          names.arg = plot_data$S,
          col = col,
          space = 1.25,
          border = NA,
          main = if(titles) "Blocking Coalitions Excess" else "",
          ylab = if(titles) "x(S) - v(S)" else "",
          las = 1,
          cex.names = 0.8,
          ylim = c(0, max(pretty(c(0, max(plot_data$excess) * 1.1)))),
          ...)
}



#### Core Points ####

#' Check Core Emptiness
#'
#' @description
#' This function determines whether the \emph{core} of a cooperative game is empty
#' by solving a Linear Programming (LP) problem. If the \emph{core} is non-empty, it returns
#' a feasible stable allocation and computes all the geometric vertices of the \emph{core} polytope.
#'
#' @details
#' The problem of finding a cost allocation \eqn{x \in \mathbb{R}^N} in the \emph{core}
#' is formulated as a Linear Programming system. The allocation must satisfy two fundamental conditions:
#'
#' \enumerate{
#' \item Efficiency: the total cost is completely allocated among the agents, i.e.,
#' \deqn{\sum_{i \in N} x_i = v(N).}
#'
#' \item Coalitional rationality: no coalition pays more than its stand-alone cost, i.e.,
#' \deqn{\sum_{i \in S} x_i \le v(S) \quad \forall S \subset N.}
#' }
#'
#' The function uses the \code{lpSolve} package to evaluate this system of equations and
#' inequalities. If a solution exists, the \emph{core} is non-empty and a feasible allocation is returned.
#'
#' Furthermore, for games with up to 6 players (\eqn{n \le 6}), the function automatically
#' utilizes the \code{rcdd} package to perform vertex enumeration. It converts the
#' system of linear inequalities into the exact geometric vertices
#' of the \emph{core}. This is particularly useful for analyzing the bounds
#' of the \emph{core} region.
#'
#' The \code{plot} method provides a visualization of the \emph{core} and the computed
#' feasible point by internally calling \code{\link{mcstCorePlot}}. It is supported
#' for games with 2, 3, or 4 players.
#'
#' @param game list; an object containing a cooperative game computed with
#' \code{method = "exact"} (e.g., from \code{\link{mcstGamePrivate}}). Default is \code{NULL}.
#' @param v numeric vector; the characteristic function values \eqn{v(S)} for all
#' \eqn{2^n - 1} coalitions, ordered by cardinality. Required if \code{game} is \code{NULL}.
#'
#' @return A list containing:
#' \itemize{
#'    \item \code{is_empty}: logical; if \code{TRUE}, the \emph{core} is empty (no stable allocation satisfies all conditions).
#'    \item \code{core_point}: numeric vector; a feasible allocation inside the \emph{core}. \code{NULL} if the \emph{core} is empty.
#'    \item \code{vertices}: a matrix where each row represents a geometric vertex of the \emph{core} polytope. \code{NULL} if the \emph{core} is empty or \eqn{n > 6}.
#'    \item \code{n}: the number of agents.
#'    \item \code{total}: the total cost of the game \eqn{v(N)}.
#'    \item \code{values}: the characteristic function values \eqn{v(S)}.
#' }
#'
#' @seealso
#' \code{\link{mcstCoreCheck}} for checking the stability of a specific allocation.
#'
#' \code{\link{mcstCorePlot}} for general \emph{core} region visualizations.
#'
#' \code{\link{mcstRules}} for an overview of the available rules and analysis tools in the package.
#'
#' @examples
#' # Example 1: Non-empty core
#' vvals1 <- c(10, 12, 15, 20, 22, 25, 30) # 3 players
#' cp1 <- mcstCorePoints(v = vvals1); cp1
#' cp1$vertices
#' plot(cp1)
#'
#' # Example 2: Empty core
#' vvals2 <- c(10, 10, 10, 15, 15, 15, 30)
#' cp2 <- mcstCorePoints(v = vvals2); cp2
#' \donttest{
#' # Example 3: Using a game object
#' C <- matrix(c(0, 12, 15, 12, 0, 4, 15, 4, 0), byrow = TRUE, ncol = 3)
#' mcstCorePoints(game = mcstGamePrivate(C))
#' }
#'
#' @concept Core Stability
#' @concept Cooperative Games
#' @concept Linear Programming
#'
#' @export

mcstCorePoints <- function(game = NULL, v = NULL) {

  if (!requireNamespace("lpSolve", quietly = TRUE))
    stop("Package 'lpSolve' required. Please install it")

  # Check inputs
  if (!is.null(game) && is.list(game)) {
    if (is.null(game$method) || game$method != "exact")
      stop("Required a 'game' object computed with method = 'exact'")

    v_game <- names(game)[grep("^v_", names(game))][1]
    vvals  <- game[[v_game]]
    m_cost <- game$total
    n      <- length(game$allocation)

  } else if (!is.null(v)) {
    vvals  <- v
    m_cost <- vvals[length(vvals)]
    n      <- log2(length(vvals) + 1)
    if (n %% 1 != 0) stop("Length of 'v' must be 2^n - 1")

  } else {
    stop("Provide either a 'game' object or 'v'")
  }


  ## Execution ##

  # Generate coalitions in cardinality order
  coals <- unlist(lapply(1:n, function(k) combn(n, k, simplify = FALSE)), recursive = FALSE)
  ncoals <- 2^n - 1

  # Build coalitional constraints: x(S) <= v(S) for all coalitions (Ax <= b)
  A <- matrix(0, nrow = ncoals, ncol = n)
  b <- numeric(ncoals)

  for (i in 1:ncoals) {
    S <- coals[[i]]
    A[i, S] <- 1
    b[i] <- vvals[i]
  }

  # Linear Programming: find any point in Core
  res <- lpSolve::lp("min", rep(0, n), rbind(rep(1, n), A),
                     c("=", rep("<=", ncoals)), c(m_cost, b))

  is_empty <- (res$status != 0)
  vertices <- NULL

  # Vertex enumeration
  if (!is_empty && n <= 6) {
    if (!requireNamespace("rcdd", quietly = TRUE)) {
      stop("Package 'rcdd' required for vertex enumeration. Please install it")
    } else {
      # H-representation matrix for rcdd. Format columns: [type, b, -A]
      # type = 1: equality (sum x = total cost)
      # type = 0: inequality (x(S) <= v(S))
      h_mat <- rbind(c(1, m_cost, -rep(1, n)), cbind(0, b, -A))

      # Convert to exact rational fractions (d2q), find vertices (scdd), and convert back to decimal (q2d)
      v_rep <- rcdd::scdd(rcdd::d2q(h_mat))
      out   <- rcdd::q2d(v_rep$output)

      # Filter valid vertices
      vertices <- out[out[, 1] == 0 & out[, 2] == 1, -c(1, 2), drop = FALSE]

      # Clean up potential floating-point duplicates
      if (!is.null(vertices)) {
        vertices <- vertices[!duplicated(round(vertices, 8)), , drop = FALSE]
      }
    }
  }
  core_point <- if (!is_empty) setNames(res$solution, as.character(seq_len(n))) else NULL


  ## Output ##

  output <- list(
    is_empty = is_empty,
    core_point = core_point,
    vertices = vertices,
    n = n,
    total = m_cost,
    values = vvals
    )

    class(output) = "mcstp_core_point"

    return(output)
}

#' @export
print.mcstp_core_point <- function(x, ...) {
  cat("-------------------------\n")
  cat(" Core Emptiness Analysis\n")
  cat("-------------------------\n")
  cat("Players (n) :", x$n, "\n")
  cat("Total Cost  :", round(x$total, 2), "\n\n")
  if (x$is_empty) {
    cat("Result: the core IS EMPTY\n")
    cat("No allocation satisfies efficiency and coalitional rationality simultaneously\n")
  } else {
    cat("Result: the core IS NOT EMPTY\n")
    cat("A feasible core allocation:\n")
    print(round(x$core_point, 2))
  }
  invisible(x)
}

#' @rdname mcstCorePoints
#' @param x the output from \code{mcstCorePoints}.
#' @param titles logical; if \code{TRUE} (default), adds a main title and subtitle/legend information to the plot.
#' @param ... additional graphical parameters passed to \code{\link{mcstCorePlot}}.
#' @export
plot.mcstp_core_point <- function(x, titles = TRUE, ...) {
  if (x$is_empty) stop("The core is empty: nothing to plot")
  if (x$n > 4) stop("Visualization is only supported for 2, 3, or 4 players")

  p <- mcstCorePlot(v = x$values, allocations = list(lpSolve = x$core_point),
                    titles = titles, ...)
  if (x$n == 4) {
    print(p)
  }

  invisible(p)
}



#### Core Plot ####

#' Plot the Core Region and Cost Allocations
#'
#' @description
#' This function generates a geometric visualization of the \emph{core} region and the
#' imputation set for cooperative games with 2, 3, or 4 players. It also allows plotting
#' multiple cost allocation rules simultaneously to visually check their stability.
#'
#' @details
#' The geometric representation of the \emph{core} and the imputation set depends on the
#' number of players \eqn{n}:
#'
#' \enumerate{
#' \item \eqn{n = 2} players: the imputation set is represented as a 1D line segment
#' where \eqn{x_1 + x_2 = v(N)}. The \emph{core} is highlighted as a sub-segment within it,
#' delimited by individual rationality.
#'
#' \item \eqn{n = 3} players: allocations are projected into 2D barycentric coordinates
#' (simplex). The imputation set forms a large triangle, and the \emph{core} is plotted inside
#' as a bounded convex polygon.
#'
#' \item \eqn{n = 4} players: the function switches to an interactive 3D plot using the
#' \code{plotly} package. The imputation set is drawn as a transparent 3D tetrahedron (mesh),
#' and the \emph{core} is represented as a solid 3D polytope inside it.
#' }
#'
#' If an allocation lies inside the green \emph{core} region, it means the allocation is
#' coalitionally stable (no group of players has incentives to defect).
#'
#' The function automatically extracts the default allocation if a \code{game} object is passed.
#' Alternatively, a custom named \code{list} of allocations can be provided to compare different
#' rules (e.g., Bird, Shapley value, folk) on the same graph.
#'
#' @param game list; an object containing a cooperative game (e.g., from \code{\link{mcstGamePrivate}}).
#' If provided and \code{allocations} is empty, its allocation and solution concept are used by default.
#' Default is \code{NULL}.
#' @param v numeric vector; the characteristic function values \eqn{v(S)} for all
#' \eqn{2^n - 1} coalitions, ordered by cardinality. Required if \code{game} is \code{NULL}.
#' @param allocations list; a named list of numeric vectors representing different cost allocations
#' to be plotted as points. Default is \code{list()}.
#' @param titles logical; if \code{TRUE} (default), adds a main title and subtitle/legend information
#' to the plot.
#'
#' @return For \eqn{n = 2} or \eqn{n = 3}, the function draws a base R plot.
#' For \eqn{n = 4}, it returns an interactive \code{plotly} object containing the 3D visualization.
#'
#' @seealso
#' \code{\link{mcstCoreCheck}} for checking stability numbers and finding blocking coalitions.
#'
#' \code{\link{mcstRules}} for an overview of the available rules and analysis tools in the package.
#'
#' @examples
#' # Example 1
#' vvals1 <- c(10, 12, 15, 20, 22, 25, 30) # 3 players: v(1), v(2), v(3), v(12), ...
#' allocs1 <- list(RuleA = c(8, 10, 12), RuleB = c(9, 10, 11), RuleC = c(7, 11, 12))
#' mcstCorePlot(v = vvals1, allocations = allocs1)
#' \donttest{
#' # Example 2 (interactive 3D plot)
#' vvals2 <- c(10, 10, 10, 10, 18, 18, 18, 18, 18,
#'             18, 25, 25, 25, 25, 32) # 4 players: v(1), v(2), v(3), v(4), v(12), ...
#' allocs2 <- list(RuleA = c(8, 8, 8, 8), RuleB = c(10, 6, 11, 5))
#' mcstCorePlot(v = vvals2, allocations = allocs2)
#' }
#'
#' @concept Core Stability
#' @concept Cooperative Games
#' @concept Visualization
#'
#' @export

mcstCorePlot <- function(game = NULL, v = NULL, allocations = list(), titles = TRUE) {

  # Check inputs
  if (!is.null(game) && is.list(game)) {
    v_game <- names(game)[grep("^v_", names(game))][1]
    vvals  <- game[[v_game]]
    m_cost <- game$total
    n <- length(game$allocation)

    if (length(allocations) == 0) {
      allocations <- list(game$allocation)
      names(allocations) <- game$sol
    }

  } else if (!is.null(v)) {
    vvals <- v
    m_cost <- vvals[length(vvals)]
    n <- log2(length(vvals) + 1)

  } else stop("Provide either a 'game' object or 'v'")

  # Calculate core vertices
  cp <- mcstCorePoints(v = vvals)
  if (cp$is_empty) stop("The core is empty: nothing to plot")

  # Colors for imputation set and core
  col_core2 <- rgb(0.45, 0.77, 0.46, 0.5); col_core3 <- "#74C476"
  col_imp2  <- rgb(0.90, 0.90, 0.90, 0.7); col_imp3  <- "#E5E5E5"

  # Colors and symbols for allocation rules
  pt_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
               "#F781BF", "#984EA3", "#A65628")
  pt_pchs <- c(15, 17, 18, 19, 16, 8, 11, 7)
  plotly_syms <- c("square", "triangle-up", "diamond", "circle",
                   "circle", "asterisk", "hexagram", "x")

  # Filter allocations to match the number of players
  allocations <- allocations[sapply(allocations, function(x) length(as.numeric(x)) == n)]


  # CASE 1: n = 2 players
  if (n == 2) {
    # Limits for the segment
    xmin <- m_cost - vvals[2]
    xmax <- vvals[1]
    all_x <- if (length(allocations) > 0) sapply(allocations, function(a) as.numeric(a)[1]) else c()
    x_range <- max(xmax, all_x) - min(xmin, all_x)
    if (x_range == 0) x_range <- 1
    x_lims <- c(min(xmin, all_x) - 0.2 * x_range, max(xmax, all_x) + 0.2 * x_range)
    y_lims <- c(m_cost - x_lims[2], m_cost - x_lims[1])

    # Plot
    old_par <- par(mar = c(4, 4, 4, 7), xpd = TRUE)
    on.exit(par(old_par))

    plot(NA, xlim = x_lims, ylim = y_lims, axes = FALSE, asp = 1,
         xlab = "", ylab = "", bty = "n")
    if (titles) {
      title(main = "Core Region & Allocations", line = 1.5)
      mtext("Imputation set (Grey) | Core (Green)", side = 3, line = 0.1, col = "darkgray", cex = 0.8)
    }

    # Imputation set
    segments(xmin, m_cost - xmin, xmax, m_cost - xmax,
             col = "grey75", lwd = 3)
    text(xmin, m_cost - xmin, labels = paste0("(", round(xmin, 2), ", ", round(m_cost - xmin, 2), ")"),
         pos = 2, cex = 0.6)
    text(xmax, m_cost - xmax, labels = paste0("(", round(xmax, 2), ", ", round(m_cost - xmax, 2), ")"),
         pos = 4, cex = 0.6)

    # Core region with vertices
    segments(xmin, m_cost - xmin, xmax, m_cost - xmax,
             col = adjustcolor(col_core2, alpha.f = 0.7), lwd = 6)
    points(c(xmin, xmax), c(m_cost - xmin, m_cost - xmax), pch = 19, cex = 0.5)

    # Individual allocation points
    if (length(allocations) > 0) {
      for (i in seq_along(allocations)) {
        points(as.numeric(allocations[[i]])[1], as.numeric(allocations[[i]])[2],
               col = pt_cols[(i - 1) %% length(pt_cols) + 1],
               pch = pt_pchs[(i - 1) %% length(pt_pchs) + 1])
        }

      legend("right", inset = c(-0.35, 0),
             legend = names(allocations),
             col = pt_cols[seq_along(allocations)],
             pch = pt_pchs[seq_along(allocations)],
             bty = "n", cex = 0.8)
    }
  }


  # CASE 2: n = 3 players
  else if (n == 3) {
    # Vertices of the imputation triangle
    imp_v <- rbind(c(vvals[1], vvals[2], m_cost - vvals[1] - vvals[2]),
                   c(vvals[1], m_cost - vvals[1] - vvals[3], vvals[3]),
                   c(m_cost - vvals[2] - vvals[3], vvals[2], vvals[3]))

    # Project 3D coordinates to 2D Barycentric coordinates
    to2d <- function(pts) {
      if (is.null(dim(pts))) pts <- matrix(pts, nrow = 1)
      cbind((pts[,1] - pts[,2]) / sqrt(2), (pts[,1] + pts[,2] - 2*pts[,3]) / sqrt(6))
    }

    imp_2d <- to2d(imp_v)
    core_2d <- to2d(cp$vertices)
    hull <- chull(core_2d[,1], core_2d[,2]) # To order core vertices

    alloc_2d <- if (length(allocations) > 0) {
      do.call(rbind, lapply(allocations, function(a) to2d(as.numeric(a))))
    } else NULL
    all_x <- c(imp_2d[,1], alloc_2d[,1])
    all_y <- c(imp_2d[,2], alloc_2d[,2])

    # Plot
    old_par <- par(mar = c(2, 2, 4, 7), xpd = TRUE)
    on.exit(par(old_par))

    plot(NULL,
         xlim = range(all_x) + c(-0.1, 0.1) * diff(range(all_x)),
         ylim = range(all_y) + c(-0.1, 0.2) * diff(range(all_y)),
         asp = 1, axes = FALSE, xlab = "", ylab = "", bty = "n")
    if (titles) {
      title(main = "Core Region & Allocations", line = 1.5)
      mtext("Imputation set (Grey) | Core (Green)", side = 3, line = 0.1, col = "darkgray", cex = 0.8)
    }

    # Imputation set
    polygon(imp_2d[,1], imp_2d[,2], border = "grey75",
            col = adjustcolor(col_imp2, alpha.f = 0.7))
    vertex_pos <- c(3, 4, 2)
    for (i in 1:3) {
      text(imp_2d[i, 1], imp_2d[i, 2],
           labels = paste0("(", round(imp_v[i, 1], 2), ", ",
                           round(imp_v[i, 2], 2), ", ",
                           round(imp_v[i, 3], 2), ")"),
           pos = vertex_pos[i], cex = 0.6)
    }

    # Core region with vertices
    polygon(core_2d[hull,1], core_2d[hull,2],
            col = adjustcolor(col_core2, alpha.f = 0.5), border = NA, lwd = 1.2)
    points(core_2d[,1], core_2d[,2], pch = 19, col = "black", cex = 0.5)

    # Individual allocation points
    if (length(allocations) > 0) {
      for (i in seq_along(allocations)) {
        p2d <- to2d(as.numeric(allocations[[i]]))
        points(p2d[1], p2d[2],
               col = pt_cols[(i - 1) %% length(pt_cols) + 1],
               pch = pt_pchs[(i - 1) %% length(pt_pchs) + 1])
      }

      legend("right", inset = c(-0.35, 0),
        legend = names(allocations),
        col = pt_cols[seq_along(allocations)],
        pch = pt_pchs[seq_along(allocations)],
        bty = "n", cex = 0.8)
    }
  }


  # CASE 3: n = 4 players
  else if (n == 4) {
    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Package 'plotly' required for 3D plots. Please install it")
    if (!requireNamespace("geometry", quietly = TRUE))
      stop("Install 'geometry'. Please install it")

    # Vertices of the imputation set
    imp_v4 <- rbind(c(vvals[1], vvals[2], vvals[3], m_cost - sum(vvals[1:3])),
                    c(vvals[1], vvals[2], m_cost - (vvals[1] + vvals[2] + vvals[4]), vvals[4]),
                    c(vvals[1], m_cost - (vvals[1] + vvals[3] + vvals[4]), vvals[3], vvals[4]),
                    c(m_cost - sum(vvals[2:4]), vvals[2], vvals[3], vvals[4]))

    # Plot
    fig <- plotly::plot_ly()

    # Imputation set
    ih <- geometry::convhulln(imp_v4[, 1:3])
    fig <- fig |> plotly::add_trace(type = "mesh3d", x = imp_v4[,1], y = imp_v4[,2], z = imp_v4[,3],
                                      i = ih[,1] - 1, j = ih[,2] - 1, k = ih[,3] - 1, intensity = 0, showscale = FALSE,
                                      colorscale = list(list(0, col_imp3), list(1, col_imp3)), opacity = 0.15, showlegend = F)

    # Core region with vertices
    vx <- cp$vertices[, 1]
    vy <- cp$vertices[, 2]
    vz <- cp$vertices[, 3]

    if (nrow(cp$vertices) >= 4) {
      ch <- geometry::convhulln(cp$vertices[, 1:3])
      fig <- fig |> plotly::add_trace(type = "mesh3d", x = vx, y = vy, z = vz,
                                      i = ch[,1] - 1, j = ch[,2] - 1, k = ch[,3] - 1, intensity = 0, showscale = FALSE,
                                      colorscale = list(list(0, col_core3), list(1, col_core3)), opacity = 0.3, showlegend = F)
    }
    fig <- fig |> plotly::add_trace(x = vx, y = vy, z = vz, type = "scatter3d", mode = "markers",
                                    marker = list(size = 1.2, color = "black", opacity = 0.6), showlegend = F)

    # Individual allocation points
    if (length(allocations) > 0) {
      for (i in seq_along(allocations)) {
        p <- as.numeric(allocations[[i]])
        fig <- fig |> plotly::add_trace(x = p[1], y = p[2], z = p[3], type = "scatter3d", mode = "markers",
                                        marker = list(size = 2.5, color = pt_cols[(i-1) %% length(pt_cols) + 1],
                                                      symbol = plotly_syms[(i-1) %% length(plotly_syms) + 1]), name = names(allocations)[i])
      }
    }

    subt <- "<span style='color:rgba(128,128,128,0.55);'>Imputation set</span> | <span style='color:rgba(49,163,84,0.55);'>Core</span>"
    fig <- fig |> plotly::layout(
      margin = list(t = 80, r = 100),
      title = list(
        text = if(titles) paste0("<b>Core Region & Allocations</b><br><span style='font-size:12px'>", subt, "</span>") else "",
        x = 0.5, y = 0.9, xanchor = "center"
      ),

      legend = list(
        y = 0.5, yanchor = "middle", x = 1.05, font = list(size = 11),
        itemwidth = 15, itemsizing = "constant", tracegroupgap = 0
      ),

      scene = list(xaxis = list(title = "P1"), yaxis = list(title = "P2"), zaxis = list(title = "P3"), aspectmode = "data")
    )
    return(fig)
  }
}



#### Stability Range ####

#' Stability Range for MCST Arcs
#'
#' @description
#' This function calculates the stability range for the arcs of a
#' minimum cost spanning tree (MCST). It determines how much the cost of an individual
#' arc can increase or decrease without changing the topology of the optimal tree.
#'
#' @details
#' The stability range of an arc \eqn{(i, j)} indicates the interval \eqn{[LB, UB]}
#' within which its original cost \eqn{c_{ij}} can vary while preserving the current minimum
#' cost spanning tree as optimal.
#'
#' At the finite boundary points of the stability interval, the minimum cost spanning tree
#' becomes non-unique due to cost ties. In these scenarios, component or path-based rules
#' (e.g., Folk, Kar) remain stable, whereas permutation-based rules such as Bird's rule or
#' Dutta-Kar's rule average allocations over all alternative optimal trees, which
#' can induce non-proportional jumps in the final cost allocation vector.
#'
#' @param C cost matrix between nodes. Accepts multiple formats (see
#' "Supported formats for \code{C}" in \code{\link{mcstBird}}).
#'
#' @return A data frame containing:
#' \itemize{
#'    \item \code{from}: the starting node of the arc.
#'    \item \code{to}: the ending node of the arc.
#'    \item \code{in_mst}: logical; if \code{TRUE}, the arc belongs to the optimal tree.
#'    \item \code{cost}: the original cost of the arc \eqn{c_{ij}}.
#'    \item \code{range}: a formatted string showing the allowed cost variation \eqn{\small [\Delta LB, \Delta UB]}.
#'    \item \code{LB}: the absolute lower bound for the arc's cost.
#'    \item \code{UB}: the absolute upper bound for the arc's cost.
#' }
#'
#' @seealso
#' \code{\link{mcstBird}}, \code{\link{mcstDuttaKar}}, \code{\link{mcstFolk}}, \code{\link{mcstKar}} for allocation rules.
#'
#' \code{\link{mcstRules}} for an overview of the available rules and analysis tools in the package.
#'
#' @examples
#' # Matrix input
#' C_mat <- matrix(c(0, 10, 15, 20,
#'                  10,  0, 25, 12,
#'                  15, 25,  0,  8,
#'                  20, 12,  8,  0), byrow = TRUE, ncol = 4)
#'
#' # Calculate stability range
#' mcstStabilityRange(C_mat)
#'
#' @concept Sensitivity Analysis
#' @concept MCSTP
#'
#' @export

mcstStabilityRange <- function(C) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)
  N0 <- rownames(C_mat)
  n0 <- length(N0)

  # Get current MST arcs and the irreducible matrix
  mst <- .get_arcs(C_mat)
  C_irred <- .get_irreducible(C_mat)

  results <- list()


  ## Execution ##

  for (k in 1:(n0 - 1)) {
    for (l in (k + 1):n0) {
      i <- N0[k]
      j <- N0[l]
      cost <- C_mat[i, j]

      # Check if the edge (i, j) is part of the current MST
      in_mst <- any((mst$i == i & mst$j == j) | (mst$i == j & mst$j == i))

      if (in_mst) { # Arc is in the MST

        C_mod <- C_mat
        C_mod[i, j] <- C_mod[j, i] <- Inf
        lb <- 0
        ub <- .get_irreducible(C_mod)[i, j]

      } else { # Arc is NOT in the MST

        lb <- C_irred[i, j]
        ub <- Inf
      }


      ## Output ##

      results[[length(results) + 1]] <- data.frame(
        from = i,
        to = j,
        in_mst = in_mst,
        cost = cost,
        range = paste0("[", round(lb - cost, 2), ", ",
                       ifelse(is.infinite(ub), "Inf", round(ub - cost, 2)), "]"),
        LB = lb,
        UB = ub,
        stringsAsFactors = FALSE
      )
    }
  }

  output <- do.call(rbind, results)
  class(output) <- c("mcstp_stability", "data.frame")

  return(output)
}

#' @export
print.mcstp_stability <- function(x, ...) {
  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}



#### Sensitivity Analysis ####

#' Sensitivity Analysis for Cost Allocations
#'
#' @description
#' This function performs a sensitivity analysis on cost allocation rules for a minimum
#' cost spanning tree (MCST) problem by varying the costs of specific arcs and observing
#' the impact on the agents' cost allocation.
#'
#' @details
#' The sensitivity analysis can be conducted under three different scenarios depending on
#' the structure of the \code{arcs} parameter and the \code{independent} flag:
#'
#' \enumerate{
#'   \item Single arc variation: when a single arc is provided, its cost varies
#'   along the designated \code{delta} interval. The output plot displays evolution lines
#'   for each agent's cost allocation relative to the absolute cost of the arc.
#'   \item Joint arc variation: when multiple arcs are provided and
#'   \code{independent = FALSE}, the same \code{delta} variation value is applied
#'   simultaneously to all selected arcs.
#'   \item Independent variation of two arcs: when exactly two arcs are provided
#'   and \code{independent = TRUE}, a grid of cross-variations is constructed using
#'   \code{delta1} and \code{delta2}. The output plot renders a selection of heatmaps
#'   (one per agent) mapping the allocations across the variation space.
#' }
#'
#' @param C cost matrix between nodes. Accepts multiple formats (see
#' "Supported formats for \code{C}" in \code{\link{mcstBird}}).
#' @param rule character string specifying the allocation rule to evaluate.
#' Must be one of \code{"bird"}, \code{"duttakar"}, \code{"folk"}, or \code{"kar"}.
#' @param arcs a vector of length 2 specifying a single arc \eqn{(i, j)}, or a list of vectors
#' specifying multiple arcs to analyze.
#' @param delta a numeric vector specifying the lower and upper limits of cost variation
#' \eqn{\small [\Delta LB, \Delta UB]}, or a list of two vectors if analyzing two independent arcs.
#' Default is \code{list(c(-5, 5), c(-3, 3))}.
#' @param step numeric value specifying the step size for the sequence of variations. Default is 1.
#' @param independent logical; if \code{TRUE}, evaluates the cross-sensitivity of exactly two arcs
#' independently, generating an evaluation grid. Default is \code{FALSE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{delta} (or \code{delta1} and \code{delta2}): the variation step values applied.
#'   \item \code{arc_cost}: the resulting absolute cost of the arc (computed only in the single arc case).
#'   \item Columns named after each agent showing the dynamic cost allocations obtained under the selected \code{rule}.
#' }
#'
#' @seealso
#' \code{\link{mcstBird}}, \code{\link{mcstDuttaKar}}, \code{\link{mcstFolk}}, \code{\link{mcstKar}} for the rule definitions.
#'
#' \code{\link{mcstStabilityRange}} for calculating static boundaries.
#'
#' \code{\link{mcstRules}} for an overview of the available rules and analysis tools in the package.
#'
#' @examples
#' # Matrix input
#' C_mat <- matrix(c(0, 10, 15, 20,
#'                   10,  0, 25, 12,
#'                   15, 25,  0,  8,
#'                   20, 12,  8,  0), byrow = TRUE, ncol = 4)
#'
#' # Example 1: Single arc sensitivity analysis
#' sens_single <- mcstSensitivity(C_mat, rule = "bird", arcs = c("0", "1"), delta = c(-5, 15))
#' plot(sens_single)
#'
#' # 2. Joint sensitivity analysis for multiple arcs
#' arcs_list <- list(c("0", "1"), c("2", "3"))
#' sens_joint <- mcstSensitivity(C_mat, rule = "folk", arcs = arcs_list, delta = c(-3, 3))
#' plot(sens_joint)
#'
#' @concept Sensitivity Analysis
#' @concept MCSTP
#'
#' @export

mcstSensitivity <- function(C, rule = c("bird", "duttakar", "folk", "kar"),
                             arcs, delta = list(c(-5, 5), c(-3, 3)), step = 1,
                             independent = FALSE) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)

  rule <- match.arg(rule)

  # Ensure delta is a list
  if (!is.list(delta)) {
    delta <- list(delta, delta)
  }

  # Generate evaluation sequences
  deltas1 <- seq(delta[[1]][1], delta[[1]][2], by = step)
  deltas2 <- seq(delta[[2]][1], delta[[2]][2], by = step)

  deltas <- deltas1
  n_deltas <- length(deltas)

  # Standardize arcs input to a list of vectors
  if (!is.list(arcs)) {
    if (is.vector(arcs) && length(arcs) == 2) {
      arcs <- list(arcs)
    } else {
      stop("Arcs must be a vector of length 2 or a list of vectors of length 2")
    }
  }


  # CASE 1: Single arc
  if (length(arcs) == 1) {
    i <- as.character(arcs[[1]][1])
    j <- as.character(arcs[[1]][2])
    orig_costs <- C_mat[i, j]

    results <- vector("list", n_deltas)

    # Evaluate the rule for each delta
    for (k in seq_along(deltas)) {
      d <- deltas[k]
      C_mod <- C_mat

      # Prevent negative costs
      C_mod[i, j] <- C_mod[j, i] <- max(0, orig_costs + d)

      # Calculate the allocation
      alloc <- switch(rule,
                      bird = mcstBird(C_mod)$e_bird,
                      duttakar = mcstDuttaKar(C_mod)$e_dk,
                      folk = mcstFolk(C_mod)$folk,
                      kar = mcstKar(C_mod)$allocation)

      results[[k]] <- c(delta = d, arc_cost = C_mod[i, j], alloc)
    }

    output <- as.data.frame(do.call(rbind, results))
    class(output) <- c("mcstp_sens", "data.frame")
    attr(output, "arc_name") <- paste0(i, "-", j)

    return(output)
  }


  # CASE 2: Multiple arcs varying simultaneously
  if (length(arcs) > 1 && !independent) {
    # Extract original costs for all selected arcs
    orig_costs <- sapply(arcs, function(e)
      C_mat[as.character(e[1]), as.character(e[2])])

    results <- vector("list", n_deltas)

    # Evaluate the rule for each delta
    for (k in seq_along(deltas)) {
      d <- deltas[k]
      C_mod <- C_mat

      # Apply the same delta to all selected arcs
      for (m in seq_along(arcs)) {
        i <- as.character(arcs[[m]][1])
        j <- as.character(arcs[[m]][2])

        # Prevent negative costs
        C_mod[i, j] <- C_mod[j, i] <- max(0, orig_costs[m] + d)
      }

      # Calculate the allocation
      alloc <- switch(rule,
                      bird = mcstBird(C_mod)$e_bird,
                      duttakar = mcstDuttaKar(C_mod)$e_dk,
                      folk = mcstFolk(C_mod)$folk,
                      kar = mcstKar(C_mod)$allocation)

      results[[k]] <- c(delta = d, alloc)
    }

    output <- as.data.frame(do.call(rbind, results))
    class(output) <- c("mcstp_sens_joint", "data.frame")
    attr(output, "arc_names") <- sapply(arcs, function(e) paste0(e[1], "-", e[2]))

    return(output)
  }


  # CASE 3: Two arcs varying independently
  if (length(arcs) > 1 && independent) {
    if (length(arcs) != 2) {
      stop("Independent sensitivity analysis is only supported for 2 arcs")
    }

    i1 <- as.character(arcs[[1]][1]); j1 <- as.character(arcs[[1]][2])
    i2 <- as.character(arcs[[2]][1]); j2 <- as.character(arcs[[2]][2])
    orig_costs1 <- C_mat[i1, j1]
    orig_costs2 <- C_mat[i2, j2]

    results <- vector("list", length(deltas1) * length(deltas2))
    idx <- 1

    # Evaluate the rule for each delta
    for (d1 in deltas1) {
      for (d2 in deltas2) {
        C_mod <- C_mat

        # Prevent negative costs
        C_mod[i1, j1] <- C_mod[j1, i1] <- max(0, orig_costs1 + d1)
        C_mod[i2, j2] <- C_mod[j2, i2] <- max(0, orig_costs2 + d2)

        # Calculate the allocation
        alloc <- switch(rule,
                        bird = mcstBird(C_mod)$e_bird,
                        duttakar = mcstDuttaKar(C_mod)$e_dk,
                        folk = mcstFolk(C_mod)$folk,
                        kar = mcstKar(C_mod)$allocation)

        results[[idx]] <- c(delta1 = d1, delta2 = d2, alloc)
        idx <- idx + 1
      }
    }

    output <- as.data.frame(do.call(rbind, results))
    class(output) <- c("mcstp_sens_indep", "data.frame")
    attr(output, "arc_names") <- c(paste0(i1, "-", j1), paste0(i2, "-", j2))

    return(output)
  }
}

#' @rdname mcstSensitivity
#' @param x the output from \code{mcstSensitivity}.
#' @param titles logical; if \code{TRUE} (default), adds informative main titles and axis labels
#' to the plots.
#' @param ... additional graphical parameters passed to the plot function.
#' @export
plot.mcstp_sens <- function(x, titles = TRUE, ...) {
  agent_cols <- 3:ncol(x)
  agent_names <- colnames(x)[agent_cols]
  n_agents <- length(agent_cols)

  colors <- rainbow(n_agents)
  shapes <- rep(c(16, 17, 15, 18, 19, 8), length.out = n_agents)
  nodes <- unlist(strsplit(attr(x, "arc_name"), "-"))

  matplot(x$arc_cost, x[, agent_cols],
          type = "l", lty = 1, col = colors,
          xlab = if(titles) bquote(c[.(nodes[1]) * .(nodes[2])] + Delta[.(nodes[1]) * .(nodes[2])]) else "",
          ylab = if(titles) "Cost Allocation" else "",
          main = if(titles) "Sensitivity Analysis" else "", ...)

  # Original state (delta = 0)
  if (any(x$delta == 0)) {
    orig_costs <- x$arc_cost[x$delta == 0][1]
    orig_alloc <- as.numeric(x[x$delta == 0, agent_cols, drop = FALSE][1, ])

    # Vertical line marking the current cost
    abline(v = orig_costs, lty = 3, col = "gray50", lwd = 0.75)

    # Adjust point sizes dynamically to handle overlaps
    sizes <- rep(1, n_agents)
    for (val in unique(orig_alloc)) {
      tied <- which(orig_alloc == val)
      if (length(tied) > 1) {
        for (k in seq_along(tied)) sizes[tied[k]] <- max(0.4, 1.5 - (k - 1) * 0.5)
      }
    }

    # Draw points in decreasing size order
    for (idx in order(sizes, decreasing = TRUE)) {
      points(orig_costs, orig_alloc[idx],
             col = colors[idx], pch = shapes[idx], cex = sizes[idx])
    }

    if (titles) {
      mtext("Current Cost", side = 3, at = orig_costs, col = "gray40",
            cex = 0.6, line = -1.5, adj = -0.15)
    }
  }

  legend("topleft", legend = paste0("P", agent_names),
         col = colors, pch = shapes,
         lty = 1, inset = c(0.025, 0),
         bty = "n", cex = 0.6, seg.len = 1)

  invisible(x)
}

#' @export
plot.mcstp_sens_joint <- function(x, titles = TRUE, ...) {
  arcs_names <- attr(x, "arc_names")
  if (is.null(arcs_names)) arcs_names <- c("ij", "kl")

  agent_cols <- 2:ncol(x)
  agent_names <- colnames(x)[agent_cols]
  n_agents <- length(agent_cols)

  colors <- rainbow(n_agents)
  shapes <- rep(c(16, 17, 15, 18, 19, 8), length.out = n_agents)
  xlab <- sapply(arcs_names, function(u) {
    nodes <- unlist(strsplit(u, "-"))
    paste0("Delta[", nodes[1], "*", nodes[2], "]")
  })

  matplot(x$delta, x[, agent_cols],
          type = "l", lty = 1, col = colors,
          xlab = if(titles) parse(text = paste(xlab, collapse = " * ', ' * ")) else "",
          ylab = if(titles) "Cost Allocation" else "",
          main = if(titles) "Sensitivity Analysis" else "", ...)

  # Original state (delta = 0)
  if (any(x$delta == 0)) {
    orig_alloc <- as.numeric(x[x$delta == 0, agent_cols, drop = FALSE][1, ])

    # Vertical line marking the current cost
    abline(v = 0, lty = 3, col = "gray50", lwd = 0.75)

    # Adjust point sizes dynamically to handle overlaps
    sizes <- rep(1, n_agents)
    for (val in unique(orig_alloc)) {
      tied <- which(orig_alloc == val)
      if (length(tied) > 1) {
        for (k in seq_along(tied)) sizes[tied[k]] <- max(0.4, 1.5 - (k - 1) * 0.5)
      }
    }

    # Draw points in decreasing size order
    for (idx in order(sizes, decreasing = TRUE)) {
      points(0, orig_alloc[idx],
             col = colors[idx], pch = shapes[idx], cex = sizes[idx])
    }

    if (titles) {
    mtext("Original State", side = 3, at = 0, col = "gray40",
          cex = 0.6, line = -1.5, adj = -0.15)
    }
  }

  legend("topleft", legend = paste0("P", agent_names),
         col = colors, pch = shapes,
         lty = 1, inset = c(0.025, 0),
         bty = "n", cex = 0.6, seg.len = 1)

  invisible(x)
}

#' @export
plot.mcstp_sens_indep <- function(x, titles = TRUE, ...) {
  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("Package 'fields' is required for independent sensitivity plotting. Please install it")
  }

  # Grid dimensions
  d1 <- unique(x$delta1)
  d2 <- unique(x$delta2)

  agent_cols <- 3:ncol(x)
  agent_names <- colnames(x)[agent_cols]
  n_agents <- length(agent_cols)
  arcs_names <- attr(x, "arc_names")

  nrows <- ceiling(sqrt(n_agents))
  ncols <- ceiling(n_agents / nrows)
  old_par <- par(mfrow = c(nrows, ncols), mar = c(3, 3, 2.5, 1.5),
                 mgp = c(2, 0.4, 0),
                 oma = c(0, 0, if (titles) 3 else 1, 0))
  on.exit(par(old_par))

  # Generate a heatmap for each agent
  for (a in seq_along(agent_cols)) {
    # Reshape the vector into a grid matrix
    mat <- matrix(x[, agent_cols[a]], nrow = length(d1), ncol = length(d2),
                  byrow = TRUE)

    nodes1 <- unlist(strsplit(arcs_names[1], "-"))
    nodes2 <- unlist(strsplit(arcs_names[2], "-"))

    # Heat Map
    fields::image.plot(d1, d2, mat, col = hcl.colors(15, "YlOrRd", rev = TRUE),
                       xlab = if(titles) bquote(Delta[.(nodes1[1]) * .(nodes1[2])]) else "",
                       ylab = if(titles) bquote(Delta[.(nodes2[1]) * .(nodes2[2])]) else "",
                       main = if(titles) paste("Agent", agent_names[a]) else "",
                       cex.axis = 0.75, tck = -0.05, cex.main = 0.85,
                       axis.args = list(cex.axis = 0.75, tck = -0.35), ...)

    # Add contour lines if there is variation in the matrix
    if (diff(range(mat)) > 0)
      contour(d1, d2, mat, add = TRUE, col = "gray20", labcex = 0.5)

    # Original state (delta = 0)
    points(0, 0, pch = 3, lwd = 1.5, cex = 1, col = "black")
  }

  if (titles) {
    mtext("Sensitivity Analysis", side = 3, outer = TRUE,
          cex = 1.2, font = 2, line = 0.25)
  }

  invisible(x)
}



#### Compare rules ####

#' Compare Allocation Rules for MCST
#'
#' @description
#' This function computes and compares the cost allocations for a minimum cost
#' spanning tree (MCST) problem using four standard rules: Bird, Dutta-Kar, folk, and Kar.
#'
#' @param C a cost matrix between nodes. Accepts multiple formats (see
#' "Supported formats for \code{C}" in \code{\link{mcstBird}}).
#'
#' @return A data frame containing the allocations
#' for each agent under the Bird, Dutta-Kar, folk, and Kar rules.
#'
#' @seealso
#' \code{\link{mcstBird}}, \code{\link{mcstDuttaKar}}, \code{\link{mcstFolk}}, \code{\link{mcstKar}} for the rule definitions.
#'
#' \code{\link{mcstRules}} for an overview of the available rules and analysis tools in the package.
#'
#' @examples
#' # Matrix input
#' C_mat <- matrix(c(0, 10, 15, 20,
#'                  10,  0, 25, 12,
#'                  15, 25,  0,  8,
#'                  20, 12,  8,  0), byrow = TRUE, ncol = 4)
#'
#' # Compare all rules
#' comparison <- mcstCompare(C_mat)
#' comparison
#'
#' # Plot the comparison
#' plot(comparison)
#'
#' @concept Comparison
#' @concept MCSTP
#'
#' @export

mcstCompare <- function(C) {

  # Convert input into a standard cost matrix C for N_0 (agents + source)
  C_mat <- .prepare_matrix(C)


  ## Execution ##

  bird <- mcstBird(C_mat)$e_bird
  dk <- mcstDuttaKar(C_mat)$e_dk
  folk <- mcstFolk(C_mat)$folk
  kar <- mcstKar(C_mat)$allocation


  ## Output ##

  output <- data.frame(
    agent = if (!is.null(names(bird))) names(bird) else seq_along(bird),
    bird = bird,
    dutta_kar = dk,
    folk = folk,
    kar = kar,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  class(output) <- c("mcstp_compare", "data.frame")

  return(output)
}

#' @export
print.mcstp_compare <- function(x, ...) {
  cat("----------------------------------\n")
  cat(" MCST Allocation Rules Comparison\n")
  cat("----------------------------------\n")
  print(as.data.frame(x), row.names = FALSE, digits = 2)
  invisible(x)
}

#' @rdname mcstCompare
#' @param x the output from \code{mcstCompare}.
#' @param titles logical; if \code{TRUE} (default), adds informative main titles and axis labels to the plot.
#' @param ... additional graphical parameters passed to \code{\link[graphics]{barplot}}.
#' @export
plot.mcstp_compare <- function(x, titles = TRUE, ...) {

  plot_mat <- t(as.matrix(x[, -1]))
  colnames(plot_mat) <- x$agent

  cols <- c("#F26D6D", "#6CB6FF", "#8AD57A", "#BE95FF")

  bp <- barplot(plot_mat,
                beside = TRUE,
                col = cols,
                space = c(0.15, 3),
                border = NA,
                ylim = c(0, max(plot_mat, na.rm = TRUE) * 1.25),
                xlab = if (titles) "Agents" else "",
                ylab = if (titles) "Cost Allocation" else "",
                main = if (titles) "Comparison of MCST Allocation Rules" else "",
                las = 1,
                ...)

  legend("topright", legend = c("Bird", "Dutta Kar", "Folk", "Kar"),
         fill = cols, bty = "n", border = NA, ncol = 4, cex = 0.7)

  invisible(x)
}

