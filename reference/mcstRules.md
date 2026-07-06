# Overview of Rules and Analysis Tools for Minimum Cost Spanning Tree Problems

This page provides an overview of the cost allocation rules implemented
in this package for minimum cost spanning tree problems (MCSTP). The
rules are organized into two main approaches: those defined directly
through algorithms and decomposition techniques, and those defined
through the analysis of associated cooperative games. It also includes a
suite of analysis tools.

## Note

When the optimal tree is not unique, rules based on Prim's algorithm
compute the symmetric average allocation over all possible permutations.
In contrast, rules based on Kruskal's and Boruvka's algorithms, as well
as the cone-wise decomposition, are unaffected by tie-breaking choices.

For game-theoretic rules, the functions in this package can compute
either the Shapley value or the Nucleolus of the associated cooperative
game.

## Rules defined through the problem

These rules are obtained directly from the network's structure and
classical MST algorithms.

*Based on Prim's algorithm:*

These rules follow the growth of a single tree starting from the source.

- [`mcstBird`](https://luciasouto13.github.io/mcstprules/reference/mcstBird.md):
  agents sequentially connect to the growing tree, paying the cost of
  the arc through which they connect to the tree.

- [`mcstDuttaKar`](https://luciasouto13.github.io/mcstprules/reference/mcstDuttaKar.md):
  a modification of Bird's rule that introduces a pivotal switch in the
  allocation cost at each step to ensure cost monotonicity.

*Based on Kruskal's algorithm:*

Based on obligation functions, these rules follow Kruskal's logic of
merging connected components using the cheapest available arcs.

- [`mcstFolk`](https://luciasouto13.github.io/mcstprules/reference/mcstFolk.md):
  the agents divide the cost of connecting their components equally (ERO
  rule).

- [`mcstOWShapley`](https://luciasouto13.github.io/mcstprules/reference/mcstOWShapley.md):
  optimistic weighted Shapley rule. A generalization of the folk rule
  that distributes costs proportionally based on asymmetric positive
  weights.

- [`mcstPWShapley`](https://luciasouto13.github.io/mcstprules/reference/mcstPWShapley.md):
  pessimistic weighted Shapley rule. An alternative weighted extension
  based on a permutation-driven obligation function.

*Based on Boruvka's algorithm:*

- [`mcstBoruvka`](https://luciasouto13.github.io/mcstprules/reference/mcstBoruvka.md):
  the agents pay the maximum possible proportion of the cheapest arc
  selected by their component. The numerical outcome is mathematically
  equivalent to the
  [`mcstFolk`](https://luciasouto13.github.io/mcstprules/reference/mcstFolk.md).

*Based on a cone-wise decomposition:*

- [`mcstCone`](https://luciasouto13.github.io/mcstprules/reference/mcstCone.md):
  a technique that decomposes any general MCSTP into a nonnegative
  combination of elementary problems (where costs are 0 or 1). It allows
  extending rules like `"folk"`, `"owshapley"`, or the Bogomolnaia and
  Moulin family (`"bogomolnaia"`) to general networks.

## Rules defined through cooperative games

These rules are defined by first associating a cooperative game (in
characteristic function form \\v\\) to the MCSTP and then applying a
solution concept (like the Shapley value or the Nucleolus).

*Cooperative games:*

- [`mcstGamePrivate`](https://luciasouto13.github.io/mcstprules/reference/mcstGamePrivate.md):
  defines the cost of a coalition as the MST cost of the subproblem
  formed by its members and the source (pessimistic approach).

- [`mcstKar`](https://luciasouto13.github.io/mcstprules/reference/mcstKar.md):
  defined as the Shapley value of the associated private game.

- [`mcstGameIrred`](https://luciasouto13.github.io/mcstprules/reference/mcstGameIrred.md):
  based on the irreducible matrix \\\bar{C}\\, where costs reflect the
  cheapest way to connect each pair of nodes through the network.

- [`mcstGameOpt`](https://luciasouto13.github.io/mcstprules/reference/mcstGameOpt.md):
  allows a coalition to connect through nodes outside the coalition at
  no cost.

- [`mcstGamePublic`](https://luciasouto13.github.io/mcstprules/reference/mcstGamePublic.md):
  assigns to each coalition the minimum cost achievable when connections
  through other agents are publicly available.

- [`mcstGameCC`](https://luciasouto13.github.io/mcstprules/reference/mcstGameCC.md):
  based on the cycle-complete matrix \\C^\*\\, which induces a concave
  cooperative game.

*Solution concepts:*

- Shapley value:

  - Based on permutations: the average of marginal contributions over
    all \\n!\\ possible orderings (\\\Pi_N\\):

    \$\$Sh_i(N, v) = \frac{1}{n!} \sum\_{\pi \in \Pi_N} \[v(Pre(i, \pi)
    \cup \\i\\) - v(Pre(i, \pi))\],\$\$

    where \\Pre(i, \pi)\\ is the set of players that precede agent \\i\\
    in the ordering \\\pi\\.

  - Based on coalitions: a weighted sum of marginal contributions over
    all possible subcoalitions \\S\\ that do not contain \\i\\:

    \$\$Sh_i(N, v) = \sum\_{S \subseteq N \setminus \\i\\} \frac{\|S\|!
    (n - \|S\| - 1)!}{n!} \[v(S \cup \\i\\) - v(S)\].\$\$

  The Shapley value is the unique solution that satisfies the properties
  of efficiency, symmetry, additivity, and the null player property.

- Nucleolus: the unique cost allocation \\x\\ that lexicographically
  minimizes the maximum excesses (complaints) of all coalitions. In a
  cost game, the excess \\e\\ of a coalition \\S\\ measures its
  dissatisfaction:

  \$\$e(S, x) = \sum\_{j \in S} x_j - v(S).\$\$

  The nucleolus ensures that the most dissatisfied coalition is as happy
  as possible. It is always efficient and lies within the Core if the
  Core is non-empty.

## Analysis tools

The package includes a comprehensive suite of tools to analyze the
stability, sensitivity, and geometric properties of the cost allocations
and the network:

- [`mcstCoreCheck`](https://luciasouto13.github.io/mcstprules/reference/mcstCoreCheck.md):
  evaluates whether a given cost allocation belongs to the core of a
  cooperative game, identifying any blocking coalitions if unstable.

- [`mcstCorePoints`](https://luciasouto13.github.io/mcstprules/reference/mcstCorePoints.md):
  determines if the core of a game is empty via Linear Programming,
  computes a feasible stable allocation, and finds the exact geometric
  vertices of the core polytope.

- [`mcstCorePlot`](https://luciasouto13.github.io/mcstprules/reference/mcstCorePlot.md):
  generates geometric visualizations of the core region and the
  imputation set for 2, 3, or 4 players, allowing for a graphical
  stability check.

- [`mcstStabilityRange`](https://luciasouto13.github.io/mcstprules/reference/mcstStabilityRange.md):
  calculates the stability range (lower and upper bounds) for the cost
  of individual arcs without altering the optimal tree topology.

- [`mcstSensitivity`](https://luciasouto13.github.io/mcstprules/reference/mcstSensitivity.md):
  performs sensitivity analysis by varying the costs of specific arcs
  (individually, jointly, or independently) and observing the dynamic
  impact on the agents' allocations.

- [`mcstCompare`](https://luciasouto13.github.io/mcstprules/reference/mcstCompare.md):
  computes and compares the cost allocations of four standard
  algorithmic rules (Bird, Dutta-Kar, folk, and Kar) simultaneously.

## References

Bergantiños, G., & Vidal-Puga, J. (2021). A review of cooperative rules
and their associated algorithms for minimum-cost spanning tree problems.
SERIEs, 12:73-100.
