# Boruvka's Rule for MCSTP

This function computes Boruvka's cost allocation (Bergantiños and Vidal-
Puga, 2011) for a minimum cost spanning tree problem \\(N_0, C)\\.
Boruvka's rule, denoted as \\\beta^{\pi}(N_0, C)\\, is defined through
Boruvka's algorithm.

## Usage

``` r
mcstBoruvka(C, draw = FALSE, titles = TRUE)
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

- `boruvka`: the Boruvka's allocation vector \\\beta^{\pi}(N_0, C)\\.

- `arcs`: a data frame of edges \\(i, j)\\ in the final tree and the
  stage at which they were added.

- `total`: the total cost of the MCST, \\m(N_0, C)\\.

- `percentage`: the share of the total cost allocated to each agent.

- `ranking`: a ranking of agents by cost (from highest to lowest; ties
  marked with \*).

- `is_unique`: logical; if `TRUE`, the MCST is unique.

## Details

The Boruvka allocation divides the cost of the optimal tree through a
step-by-step process over \\\gamma\\ stages. At each stage \\s\\, agents
pay the maximum possible proportion \\p^s\\ of the cheapest arc
associated with their connected component.

Let \\A^s\\ be the set of non-completely paid arcs, \\\varrho\_{ij}^s\\
the proportion of the cost of arc \\(i,j)\\ already paid, and
\\N\_{ij}^s\\ the set of agents paying for arc \\(i,j)\\ in stage \\s\\.
The proportion \\p^s\\ is given by:

\$\$p^s = \min \left\\ \dfrac{1 - \varrho\_{ij}^{s-1}}{\|N\_{ij}^s\|} :
(i,j) \in A^{s-1}, ~ N\_{ij}^s \neq \emptyset \right\\.\$\$

Let \\a_i^s\\ be the arc partially paid by agent \\i\\ in stage \\s\\.
The cost paid by agent \\i\\ in this stage is:

\$\$f_i^s = p^s c\_{a_i^s}.\$\$

The process finishes in \\\gamma\\ stages when the tree is completely
paid, i.e., \\\sum\_{s=1}^\gamma p^s = 1\\. The Boruvka allocation for
each agent \\i \in N\\ is given by:

\$\$\beta_i^{\pi}(N_0, C) = \displaystyle{\sum\_{s=1}^\gamma f_i^s.}\$\$

## Note

Boruvka's allocation is uniquely determined even when the minimum cost
spanning tree is not unique. Since the rule depends on the evolution of
the connected components rather than the specific arcs selected, any
tie-breaking selection in Boruvka's algorithm yields the same cost
distribution.

Furthermore, as demonstrated by Bergantiños and Vidal-Puga (2011), for
any order \\\pi\\ used in the algorithm, Boruvka's allocation is
equivalent to the
[`mcstFolk`](https://luciasouto13.github.io/mcstprules/reference/mcstFolk.md).

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

Bergantiños G, Vidal-Puga J (2011) The folk solution and Boruvka’s
algorithm in minimum cost spanning tree problems. Discrete Appl Math
159(12):1279–1283

Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules and
their associated algorithms for minimum-cost spanning tree problems.
SERIEs, 12:73-100

## See also

[`mcstFolk`](https://luciasouto13.github.io/mcstprules/reference/mcstFolk.md)
for an equivalent rule based on Kruskal's algorithm.

[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md)
for an overview of the available rules and analysis tools in the
package.

## Examples

``` r
# Simple vector input
mcstBoruvka(c(12, 15, 20, 4, 6, 8), draw = TRUE)

#> 1 2 3 
#> 7 7 8 

# Input with infinite costs (disconnected nodes)
C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
           Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
mcstBoruvka(C_inf)
#> 1 2 3 4 5 6 
#> 5 7 6 6 6 9 

# Matrix input
C_mat <- matrix(c(0, 12, 15, 12,
                 12,  0,  4,  6,
                 15,  4,  0,  8,
                 12,  6,  8,  0), byrow = TRUE, ncol = 4)
mcstBoruvka(C_mat, draw = TRUE, titles = FALSE)

#> 1 2 3 
#> 7 7 8 
```
