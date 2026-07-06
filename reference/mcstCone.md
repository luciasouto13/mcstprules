# Cone-wise Decomposition Rules for MCSTP

This function computes cost allocations for minimum cost spanning tree
problems \\(N_0, C)\\ using the cone-wise decomposition introduced by
Norde et al. (2004). This technique allows extending rules defined for
elementary MCSTPs to any MCSTP.

## Usage

``` r
mcstCone(
  C,
  rule = c("folk", "owshapley", "bogomolnaia"),
  weights = NULL,
  lambda = 1
)
```

## Arguments

- C:

  cost matrix between nodes. Accepts multiple formats (see "Supported
  formats for `C`" below).

- rule:

  a character string indicating the rule \\R\\ to be applied to the
  elementary MCSTPs. One of `"folk"` (default), `"owshapley"`, or
  `"bogomolnaia"`.

- weights:

  a numeric vector of strictly positive weights for each agent in \\N\\.
  The length of the vector must match the number of agents; only if
  `rule = "owshapley"`.

- lambda:

  a non-negative numeric parameter (`1` by default); only if
  `rule = "bogomolnaia"`.

## Value

A list containing:

- `allocations`: the cost allocation vector \\R(N_0, C)\\ for the chosen
  rule \\R\\.

- `rule`: the rule \\R\\ applied for the allocation.

- `total`: the total cost of the MCST, \\m(N_0, C)\\.

- `percentage`: the share of the total cost allocated to each agent.

- `ranking`: a ranking of agents by cost (from highest to lowest; ties
  marked with \*).

- `decomposition`: a matrix showing the allocation at each cost level of
  the decomposition.

- `weights`: the weight vector used for the allocation; only if
  `rule = "owshapley"`.

- `lambda`: the lambda value used for the allocation; only if
  `rule = "bogomolnaia"`.

## Details

Any MCSTP can be written as a nonnegative combination of elementary MCST
\\C = \sum\_{q=1}^{m(C)} x^q C^q\\, where the costs of the arcs in
\\C^q\\ are 0 or 1. A rule \\R\\ is extended as:

\$\$R(N_0, C) = \displaystyle{\sum\_{q=1}^{m(C)} x^q R(N_0, C^q).}\$\$

Let \\g^q\\ be the network at stage \\q \in \\1, \dots, m(C)\\\\, and
\\S_i^q := S(P(g^q), i)\\ the connected component containing agent \\i\\
in the graph induced by zero-cost arcs. For each elementary problem
\\(N_0, C^q)\\, the agents in a component \\S\\ not containing the
source (0) divide the cost of connecting to the source according to the
chosen rule:

- Folk rule (`"folk"`): \$\$F(N_0, C^q) = \begin{cases}
  \dfrac{1}{\|S_i^q\|} & \text{if } 0 \notin S_i^q \\ 0 &
  \text{otherwise}. \end{cases}\$\$

- Optimistic weighted Shapley rule (`"owshapley"`):
  \$\$f^{\varrho^{ow}}(N_0, C^q) = \begin{cases} \dfrac{w_i}{\sum\_{j
  \in S_i^q} w_j} & \text{if } 0 \notin S_i^q \\ 0 & \text{otherwise}.
  \end{cases}\$\$

- Bogomolnaia and Moulin family (`"bogomolnaia"`):

  For \\\lambda \in \[0, +\infty)\\, \$\$R^\lambda(N_0, C^q) =
  \begin{cases} \dfrac{\lambda^{\delta_i}}{\sum\_{j \in S_i^q}
  \lambda^{\delta_j}} & \text{if } 0 \notin S_i^q \\ 0 &
  \text{otherwise}, \end{cases}\$\$ where \\\delta_i\\ denotes the
  number of non-null arcs in \\C^q\\ containing agent \\i\\.

  For \\\lambda = +\infty\\, \$\$R^\lambda(N_0, C^q) = \begin{cases}
  \arg\max\_{j \in S_i^q} \delta_j & \text{if } 0 \notin S_i^q \\ 0 &
  \text{otherwise}, \end{cases}\$\$ i.e., agents attaining the maximum
  value of \\\delta_j\\ share the cost equally.

## Note

This function focuses on the analytical and numerical decomposition of
the cost allocation. For graphical representations, use the direct rule
functions (e.g.,
[`mcstFolk`](https://luciasouto13.github.io/mcstprules/reference/mcstFolk.md)).

The allocation for `rule = "folk"` is also obtained when all weights are
equal in `rule = "owshapley"`, or when \\\lambda = 1\\ in
`rule = "bogomolnaia"`.

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

Norde H, Moretti S, Tijs S (2004) Minimum cost spanning tree games and
population monotonic allocation schemes. Eur J Oper Res 154(1):84–97

## See also

[`mcstFolk`](https://luciasouto13.github.io/mcstprules/reference/mcstFolk.md),
[`mcstOWShapley`](https://luciasouto13.github.io/mcstprules/reference/mcstOWShapley.md)
for direct implementations of these rules.

[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md)
for an overview of the available rules and analysis tools in the
package.

## Examples

``` r
C <- c(12, 15, 20, 4, 6, 8)

# Folk rule
mcstCone(C, rule = "folk")
#> 1 2 3 
#> 7 7 8 

# Optimistic weighted Shapley rule
mcstCone(C, rule = "owshapley", weights = c(2, 1, 1))
#> Weights: 2 1 1
#>    1    2    3 
#> 8.33 6.17 7.50 

# Bogomolnaia and Moulin family
mcstCone(C, rule = "bogomolnaia", lambda = 2)
#> Lambda: 2
#>    1    2    3 
#> 6.73 7.13 8.13 
mcstCone(C, rule = "bogomolnaia", lambda = Inf)
#> Lambda: Inf
#>    1    2    3 
#> 6.33 7.33 8.33 
```
