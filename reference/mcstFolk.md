# Folk Rule for MCSTP

This function computes folk cost allocation (Feltkamp et al., 1994) for
a minimum cost spanning tree problem \\(N_0, C)\\. The folk rule,
denoted as \\F(N_0, C)\\, also known as the ERO (Equal Remaining
Obligation) rule, is defined through Kruskal's algorithm using a
function with specific properties, called the obligation function
\\o_i(S)\\.

## Usage

``` r
mcstFolk(C, draw = FALSE, titles = TRUE)
```

## Arguments

- C:

  cost matrix between nodes. Accepts multiple formats (see "Supported
  formats for `C`" below).

- draw:

  logical; if `TRUE`, plots the network highlighting an optimal tree in
  red, indicating the stage at which each arc is added (in brackets) and
  the cost allocated to each agent (in parentheses below the node).

- titles:

  logical; if `TRUE` (default), adds a main title specifying the
  allocation rule and a subtitle with the algorithm used and the total
  network cost.

## Value

A list containing:

- `folk`: the folk allocation vector \\F(N_0, C)\\.

- `arcs`: a data frame of edges \\(i, j)\\ in the final tree and the
  stage at which they were added.

- `total`: the total cost of the MCST, \\m(N_0, C)\\.

- `percentage`: the share of the total cost allocated to each agent.

- `ranking`: a ranking of agents by cost (from highest to lowest; ties
  marked with \*).

- `is_unique`: logical; if `TRUE`, the MCST is unique.

## Details

Following Kruskal's algorithm, let \\g^p\\ be the network at stage \\p
\in \\1, \dots, n\\\\, and \\S_i^p := S(P(g^p), i)\\ the connected
component containing agent \\i\\. The folk rule is based on the
obligation function:

\$\$o_i(S) = \dfrac{1}{\|S\|}.\$\$

At each stage \\p\\, let \\c\_{i^p j^p}\\ be the cost of the arc added
to the network. The folk allocation for each agent \\i \in N\\ is given
by:

\$\$F_i(N_0, C) = \displaystyle{ \sum\_{p=1}^n c\_{i^p j^p} \left(
o_i(S_i^{p-1}) - o_i(S_i^p) \right)}.\$\$

## Note

The folk allocation is uniquely determined even when the minimum cost
spanning tree is not unique. Since the rule depends on the evolution of
the connected components rather than the specific arcs selected, any
tie-breaking selection in Kruskal's algorithm yields the same cost
distribution.

## Supported formats for `C`

Note: In all cases, `Inf` is accepted for disconnected node pairs. The
first row/column is assumed to be the source node (0).

The function accepts the following formats for the cost input:

- Numeric vector:

  The lower triangle of a symmetric cost matrix, excluding the diagonal.
  Length must be \\n(n+1)/2\\, where \\n\\ is the number of agents
  (excluding the source).

- Square matrix / Data frame:

  A standard adjacency matrix where `C[i, j]` represents the cost
  between node \\i\\ and node \\j\\.

- Edge list (matrix or data frame):

  An alternative to the full matrix. It must contain columns named
  `"from"` and `"to"`, plus a third column for the weights/costs.

- 'igraph' object:

  A weighted undirected graph. See the
  [`igraph`](https://r.igraph.org/reference/aaa-igraph-package.html)
  package for details on creating these objects.

- File path:

  A string pointing to a `.csv`, `.xlsx`, or `.xls` file. The file can
  contain either a square matrix (without headers) or an edge list (with
  `"from"` and `"to"` headers).

## References

Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules and
their associated algorithms for minimum-cost spanning tree problems.
SERIEs, 12:73-100

Feltkamp V, Tijs S, Muto S (1994) On the irreducible core and the equal
remaining obligations rule of minimum cost spanning extension problems.
Technical Report 106, CentER DP 1994, Tilburg University, The
Netherlands

## See also

[`mcstOWShapley`](https://luciasouto13.github.io/mcstprules/reference/mcstOWShapley.md),
[`mcstPWShapley`](https://luciasouto13.github.io/mcstprules/reference/mcstPWShapley.md)
for other rules based on Kruskal's algorithm.

[`mcstBoruvka`](https://luciasouto13.github.io/mcstprules/reference/mcstBoruvka.md)
for an equivalent rule based on Boruvka's algorithm.

[`mcstCone`](https://luciasouto13.github.io/mcstprules/reference/mcstCone.md)
with `rule = "folk"` for the equivalent implementation based on
cone-wise decomposition.

[`mcstGameIrred`](https://luciasouto13.github.io/mcstprules/reference/mcstGameIrred.md),
[`mcstGameOpt`](https://luciasouto13.github.io/mcstprules/reference/mcstGameOpt.md)
for cooperative games whose Shapley value coincides with the folk
solution.

[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md)
for an overview of the available rules and analysis tools in the
package.

## Examples

``` r
# Simple vector input
mcstFolk(c(12, 15, 20, 4, 6, 8), draw = TRUE)

#> 1 2 3 
#> 7 7 8 

# Input with infinite costs (disconnected nodes)
C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
           Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
mcstFolk(C_inf)
#> 1 2 3 4 5 6 
#> 5 7 6 6 6 9 

# Matrix input
C_mat <- matrix(c(0, 12, 15, 12,
                 12,  0,  4,  6,
                 15,  4,  0,  8,
                 12,  6,  8,  0), byrow = TRUE, ncol = 4)
mcstFolk(C_mat, draw = TRUE, titles = FALSE)

#> 1 2 3 
#> 7 7 8 
```
