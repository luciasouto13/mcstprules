# Cost Allocations based on the Optimistic Game for MCSTP

This function computes the cost allocation for the optimistic game
(Bergantiños and Vidal-Puga, 2007b) associated with a minimum cost
spanning tree problem \\(N_0, C)\\, using a specified solution concept.
The optimistic game, denoted as \\(N, v_C^o)\\ \\(\\or simply \\(N,
v^o)\\\\)\\, evaluates the cost of each coalition \\S \subseteq N\\
assuming that agents in \\N \setminus S\\ are already connected. Thus,
agents in \\S\\ can connect to the source through agents in \\N
\setminus S\\ for free.

## Usage

``` r
mcstGameOpt(
  C,
  sol = c("shapley", "nucleolus"),
  draw = FALSE,
  which.plot = c("main", "details"),
  titles = TRUE,
  method = c("exact", "montecarlo"),
  nsim = 10000
)
```

## Arguments

- C:

  cost matrix between nodes. Accepts multiple formats (see "Supported
  formats for `C`" below).

- sol:

  character; the solution concept to apply. Options are `"shapley"`
  (default), or `"nucleolus"` (only for `method = "exact"`).

- draw:

  logical; if `TRUE`, plots the network highlighting an optimal tree in
  red (constructed using Prim's algorithm), indicating the cost
  allocated to each agent in brackets below the node. For \\n \le 3\\
  and `method = "exact"`, it also displays the optimal network for each
  subcoalition \\S \subset N\\.

- which.plot:

  character string indicating which plots to display. If `"details"`
  (only for \\n \le 3\\ and `method = "exact"`), displays the optimal
  network for each subcoalition \\S \subset N\\; if `"main"`, the
  allocation for \\N\\ is plotted. Default is `c("main", "details")`.

- titles:

  logical; if `TRUE` (default), adds a main title specifying the
  cooperative game and a subtitle with the solution concept, method
  (exact or Monte Carlo) and the total network cost.

- method:

  character string specifying the calculation method. If `"exact"`
  (default), computes the exact allocation (Shapley value or nucleolus)
  and the full characteristic function; if `"montecarlo"`, estimates an
  approximation of the Shapley value through random sampling using
  `nsim` permutations. It uses lazy evaluation to avoid computing all
  possible coalitions (recommended for large networks).

- nsim:

  integer; the number of permutations to sample when using
  `method = "montecarlo"`. Default is 10,000.

## Value

A list containing:

- `coalitions`: a character vector representing all possible coalitions
  \\S \subseteq N\\, ordered by cardinality. If `method = "montecarlo"`,
  returns a string indicating that lazy evaluation was used.

- `C_opt`: a list of the optimistic cost matrices \\C^T\\ for each
  coalition \\S \subseteq N\\. If `method = "montecarlo"`, returns a
  string indicating that lazy evaluation was used.

- `v_o`: a numeric vector containing the characteristic function values
  \\v^o(S)\\ for all coalitions. If `method = "montecarlo"`, returns a
  string indicating that lazy evaluation was used.

- `contributions`: a data frame with the marginal contributions for the
  sampled permutations. If `method = "exact"`, returns a string
  indicating that this table is not used.

- `allocation`: a numeric vector containing the computed cost allocation
  for each agent.

- `total`: the total cost of the MCST, \\v^o(N)\\.

- `nperms`: the number of permutations `nsim` used in the calculation.
  If `method = "exact"`, returns a string indicating that it was not
  computed.

- `percentage`: the share of the total cost allocated to each agent.

- `ranking`: ranking of agents by cost (from highest to lowest; ties
  marked with \*).

- `method`: the method used for the calculation (`"exact"` or
  `"montecarlo"`).

- `sol`: the solution concept used (`"shapley"` or `"nucleolus"`).

## Details

Let \\T := N \setminus S\\. The optimistic game \\(N, v_C^o)\\ is a
cooperative game with transferable utility (TU game) defined for each
coalition \\S \subseteq N\\ as:

\$\$v_C^o(S) = m(S_0, C^T),\$\$

where \\m(S_0, C^T)\\ is the cost of the minimum spanning tree of the
subproblem induced by the agents in \\S\\ and the source, using the
modified cost matrix \\C^T\\. In this matrix, the cost between agents in
\\S\\ remains the same, i.e., \\c\_{ij}^T = c\_{ij}\\ for all \\i,j \in
S\\, but the cost to connect any agent \\i \in S\\ to the source is
updated as:

\$\$c\_{0i}^T = \min\_{j \in T \cup \\0\\} c\_{ji}.\$\$

Under `method = "exact"`, the function builds the full characteristic
function (all \\2^n - 1\\ values). Depending on the chosen solution
concept, it computes the allocation via
[`shapleyValue`](https://rdrr.io/pkg/CoopGame/man/shapleyValue.html)
(which relies on the coalition-based formula) or
[`nucleolus`](https://rdrr.io/pkg/CoopGame/man/nucleolus.html).

For larger networks, `method = "montecarlo"` provides an approximation
of the Shapley value through random sampling using \\n\_{sim}\\
permutations. This approach implements lazy evaluation based on the
permutation formula of the Shapley value: instead of computing the MST
cost for all \\2^n - 1\\ possible subsets, it only evaluates the
specific subcoalitions generated by the random permutations.

This is highly useful when the number of agents is large, as it
significantly reduces memory usage and computation time while providing
a very accurate approximation of the exact allocation. For networks with
few players, the exact method is preferred to avoid potential
overestimations.

For more details on these formulas, see
[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md).

## Note

The optimistic game associated with any MCST problem coincides with the
optimistic game associated with its irreducible form, i.e.,

\$\$v_C^o = v\_{\bar{C}}^o.\$\$

If \\(N_0, \bar{C})\\ is irreducible, the private and optimistic games,
\\v^p\\ and \\v^o\\, are dual: \\v^p(S) + v^o(N \setminus S) = m(N_0,
C)\\ for all \\S \subset N\\.

Furthermore, \\Sh(N, v_C^o) = Sh(N, v_C^i)\\, so applying the Shapley
value to \\v^o\\ provides an alternative way of obtaining the
[`mcstFolk`](https://luciasouto13.github.io/mcstprules/reference/mcstFolk.md)
(Bergantiños and Vidal-Puga, 2007b).

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

Bergantiños G, Vidal-Puga J (2007b) The optimistic TU game in minimum
cost spanning tree problems. Int J Game Theory 36(2):223–239

Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules and
their associated algorithms for minimum-cost spanning tree problems.
SERIEs, 12:73-100

## See also

[`mcstGamePrivate`](https://luciasouto13.github.io/mcstprules/reference/mcstGamePrivate.md)
for a game based on a pessimistic approach.

[`mcstGameIrred`](https://luciasouto13.github.io/mcstprules/reference/mcstGameIrred.md)
for an equivalent game under the Shapley value.

[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md)
for an overview of the available rules and analysis tools in the
package.

## Examples

``` r
# Simple vector input
mcstGameOpt(c(12, 15, 20, 4, 6, 8), sol = "shapley", draw = TRUE)


#>       S v^o(S)
#> 1   {1}      4
#> 2   {2}      4
#> 3   {3}      6
#> 4 {1,2}     10
#> 5 {1,3}     10
#> 6 {2,3}     10
#> 7     N     22
#> 
#> Solution: Shapley 
#> 1 2 3 
#> 7 7 8 

# Input with infinite costs (disconnected nodes)
C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
           Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
mcstGameOpt(C_inf, sol = "shapley")
#>              S v^o(S)
#> 1          {1}      5
#> 2          {2}      7
#> 3          {3}      5
#> 4          {4}      5
#> 5          {5}      6
#> 6          {6}      9
#> 7        {1,2}     12
#> 8        {1,3}     10
#> 9        {1,4}     10
#> 10       {1,5}     11
#> 11       {1,6}     14
#> 12       {2,3}     12
#> 13       {2,4}     12
#> 14       {2,5}     13
#> 15       {2,6}     16
#> 16       {3,4}     12
#> 17       {3,5}     11
#> 18       {3,6}     14
#> 19       {4,5}     11
#> 20       {4,6}     14
#> 21       {5,6}     15
#> 22     {1,2,3}     17
#> 23     {1,2,4}     17
#> 24     {1,2,5}     18
#> 25     {1,2,6}     21
#> 26     {1,3,4}     17
#> 27     {1,3,5}     16
#> 28     {1,3,6}     19
#> 29     {1,4,5}     16
#> 30     {1,4,6}     19
#> 31     {1,5,6}     20
#> 32     {2,3,4}     19
#> 33     {2,3,5}     18
#> 34     {2,3,6}     21
#> 35     {2,4,5}     18
#> 36     {2,4,6}     21
#> 37     {2,5,6}     22
#> 38     {3,4,5}     18
#> 39     {3,4,6}     21
#> 40     {3,5,6}     20
#> 41     {4,5,6}     20
#> 42   {1,2,3,4}     24
#> 43   {1,2,3,5}     23
#> 44   {1,2,3,6}     26
#> 45   {1,2,4,5}     23
#> 46   {1,2,4,6}     26
#> 47   {1,2,5,6}     27
#> 48   {1,3,4,5}     23
#> 49   {1,3,4,6}     26
#> 50   {1,3,5,6}     25
#> 51   {1,4,5,6}     25
#> 52   {2,3,4,5}     25
#> 53   {2,3,4,6}     28
#> 54   {2,3,5,6}     27
#> 55   {2,4,5,6}     27
#> 56   {3,4,5,6}     27
#> 57 {1,2,3,4,5}     30
#> 58 {1,2,3,4,6}     33
#> 59 {1,2,3,5,6}     32
#> 60 {1,2,4,5,6}     32
#> 61 {1,3,4,5,6}     32
#> 62 {2,3,4,5,6}     34
#> 63           N     39
#> 
#> Solution: Shapley 
#> 1 2 3 4 5 6 
#> 5 7 6 6 6 9 

# Matrix input
C_mat <- matrix(c(0, 12, 15, 12,
                 12,  0,  4,  6,
                 15,  4,  0,  8,
                 12,  6,  8,  0), byrow = TRUE, ncol = 4)
mcstGameOpt(C_mat, draw = TRUE, which.plot = "main", titles = FALSE)

#>       S v^o(S)
#> 1   {1}      4
#> 2   {2}      4
#> 3   {3}      6
#> 4 {1,2}     10
#> 5 {1,3}     10
#> 6 {2,3}     10
#> 7     N     22
#> 
#> Solution: Shapley 
#> 1 2 3 
#> 7 7 8 
```
