#' Overview of Rules and Analysis Tools for Minimum Cost Spanning Tree Problems
#'
#' @description
#' This page provides an overview of the cost allocation rules implemented in this
#' package for minimum cost spanning tree problems (MCSTP).
#' The rules are organized into two main approaches: those defined directly through
#' algorithms and decomposition techniques, and those defined through the analysis
#' of associated cooperative games. It also includes a suite of analysis tools.
#'
#' @section Rules defined through the problem:
#'
#' These rules are obtained directly from the network's structure and
#' classical MST algorithms.
#'
#'
#' \emph{Based on Prim's algorithm:}
#'
#' These rules follow the growth of a single tree starting from the source.
#' \itemize{
#'   \item \code{\link{mcstBird}}: agents sequentially connect to the growing
#'   tree, paying the cost of the arc through which they connect to the tree.
#'   \item \code{\link{mcstDuttaKar}}: a modification of Bird's rule that introduces
#'   a pivotal switch in the allocation cost at each step to ensure cost monotonicity.
#' }
#'
#'
#' \emph{Based on Kruskal's algorithm:}
#'
#' Based on obligation functions, these rules follow Kruskal's logic of
#' merging connected components using the cheapest available arcs.
#' \itemize{
#'   \item \code{\link{mcstFolk}}: the agents divide the cost of connecting
#'   their components equally (ERO rule).
#'   \item \code{\link{mcstOWShapley}}: optimistic weighted Shapley rule. A
#'   generalization of the folk rule that distributes costs proportionally
#'   based on asymmetric positive weights.
#'   \item \code{\link{mcstPWShapley}}: pessimistic weighted Shapley rule. An
#'   alternative weighted extension based on a permutation-driven obligation function.
#' }
#'
#'
#' \emph{Based on Boruvka's algorithm:}
#'
#' \itemize{
#'   \item \code{\link{mcstBoruvka}}: the agents pay the maximum possible proportion
#'   of the cheapest arc selected by their component. The numerical outcome
#'   is mathematically equivalent to the \code{\link{mcstFolk}}.
#' }
#'
#'
#' \emph{Based on a cone-wise decomposition:}
#'
#' \itemize{
#'   \item \code{\link{mcstCone}}: a technique that decomposes any general
#'   MCSTP into a nonnegative combination of elementary problems (where costs
#'   are 0 or 1). It allows extending rules like \code{"folk"}, \code{"owshapley"},
#'   or the Bogomolnaia and Moulin family (\code{"bogomolnaia"}) to general networks.
#' }
#'
#'
#' @section Rules defined through cooperative games:
#'
#' These rules are defined by first associating a cooperative game (in
#' characteristic function form \eqn{v}) to the MCSTP and then applying a solution
#' concept (like the Shapley value or the Nucleolus).
#'
#'
#' \emph{Cooperative games:}
#'
#' \itemize{
#'   \item \code{\link{mcstGamePrivate}}: defines the cost of a coalition as the
#'   MST cost of the subproblem formed by its members and the source (pessimistic
#'   approach).
#'   \item \code{\link{mcstKar}}: defined as the Shapley value of the associated private game.
#'   \item \code{\link{mcstGameIrred}}: based on the irreducible matrix \eqn{\bar{C}},
#'   where costs reflect the cheapest way to connect each pair of nodes through
#'   the network.
#'   \item \code{\link{mcstGameOpt}}: allows a coalition to connect through nodes
#'   outside the coalition at no cost.
#'   \item \code{\link{mcstGamePublic}}: assigns to each coalition the minimum cost
#'   achievable when connections through other agents are publicly available.
#'   \item \code{\link{mcstGameCC}}: based on the cycle-complete matrix \eqn{C^*},
#'   which induces a concave cooperative game.
#' }
#'
#'
#' \emph{Solution concepts:}
#'
#' \itemize{
#'   \item Shapley value:
#'
#'   \itemize{
#'   \item Based on permutations: the average of marginal contributions
#'   over all \eqn{n!} possible orderings (\eqn{\Pi_N}):
#'
#'   \deqn{Sh_i(N, v) = \frac{1}{n!} \sum_{\pi \in \Pi_N} [v(Pre(i, \pi) \cup \{i\}) - v(Pre(i, \pi))],}
#'
#'   where \eqn{Pre(i, \pi)} is the set of players that precede agent \eqn{i}
#'   in the ordering \eqn{\pi}.
#'
#'   \item Based on coalitions: a weighted sum of marginal contributions
#'   over all possible subcoalitions \eqn{S} that do not contain \eqn{i}:
#'
#'   \deqn{Sh_i(N, v) = \sum_{S \subseteq N \setminus \{i\}} \frac{|S|! (n - |S| - 1)!}{n!} [v(S \cup \{i\}) - v(S)].}
#'   }
#'
#'   The Shapley value is the unique solution that satisfies the properties of
#'   efficiency, symmetry, additivity, and the null player property.
#'
#'
#'   \item Nucleolus: the unique cost allocation \eqn{x} that lexicographically minimizes the
#'   maximum excesses (complaints) of all coalitions. In a cost game, the excess \eqn{e} of a coalition
#'   \eqn{S} measures its dissatisfaction:
#'
#'   \deqn{e(S, x) = \sum_{j \in S} x_j - v(S).}
#'
#'   The nucleolus ensures that the most dissatisfied coalition is as happy as
#'   possible. It is always efficient and lies within the Core if the Core is
#'   non-empty.
#' }
#'
#'
#' @section Analysis tools:
#'
#' The package includes a comprehensive suite of tools to analyze the stability,
#' sensitivity, and geometric properties of the cost allocations and the network:
#'
#' \itemize{
#'   \item \code{\link{mcstCoreCheck}}: evaluates whether a given cost allocation belongs
#'   to the core of a cooperative game, identifying any blocking coalitions if unstable.
#'   \item \code{\link{mcstCorePoints}}: determines if the core of a game is empty via
#'   Linear Programming, computes a feasible stable allocation, and finds the exact geometric
#'   vertices of the core polytope.
#'   \item \code{\link{mcstCorePlot}}: generates geometric visualizations of the core
#'   region and the imputation set for 2, 3, or 4 players, allowing for a graphical stability check.
#'   \item \code{\link{mcstStabilityRange}}: calculates the stability range (lower and
#'   upper bounds) for the cost of individual arcs without altering the optimal tree topology.
#'   \item \code{\link{mcstSensitivity}}: performs sensitivity analysis by varying the costs
#'   of specific arcs (individually, jointly, or independently) and observing the dynamic impact
#'   on the agents' allocations.
#'   \item \code{\link{mcstCompare}}: computes and compares the cost allocations of four
#'   standard algorithmic rules (Bird, Dutta-Kar, folk, and Kar) simultaneously.
#' }
#'
#' @note
#' When the optimal tree is not unique, rules based on Prim's algorithm compute
#' the symmetric average allocation over all possible permutations.
#' In contrast, rules based on Kruskal's and Boruvka's algorithms, as well as the
#' cone-wise decomposition, are unaffected by tie-breaking choices.
#'
#' For game-theoretic rules, the functions in this package can compute either
#' the Shapley value or the Nucleolus of the associated cooperative game.
#'
#' @references
#' Bergantiños, G., & Vidal-Puga, J. (2021). A review of cooperative rules
#' and their associated algorithms for minimum-cost spanning tree problems.
#' SERIEs, 12:73-100.
#'
#' @concept Algorithmic Rules
#' @concept Cooperative Games
#' @concept Analysis Tools
#' @concept MCSTP
#'
#' @name mcstRules
NULL

