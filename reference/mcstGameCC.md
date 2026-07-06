# Cost Allocations based on the Cycle-complete Game for MCSTP

This function computes the cost allocation for the cycle-complete game
(Trudeau, 2012) associated with a minimum cost spanning tree problem
\\(N_0, C)\\, using a specified solution concept. The cycle-complete
game, denoted as \\(N, v_C^c)\\ \\(\\or simply \\(N, v^c)\\\\)\\, is
defined as the private game associated with the cycle-complete network
\\(N_0, C^\*)\\.

## Usage

``` r
mcstGameCC(
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

- `C_cc`: the cycle-complete cost matrix \\C^\*\\.

- `v_c`: a numeric vector containing the characteristic function values
  \\v^c(S)\\ for all coalitions. If `method = "montecarlo"`, returns a
  string indicating that lazy evaluation was used.

- `contributions`: a data frame with the marginal contributions for the
  sampled permutations. If `method = "exact"`, returns a string
  indicating that this table is not used.

- `allocation`: a numeric vector containing the computed cost allocation
  for each agent.

- `total`: the total cost of the MCST, \\v^c(N)\\.

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

The cycle-complete game \\(N, v_C^c)\\ is a cooperative game with
transferable utility (TU game) defined for each coalition \\S \subseteq
N\\ as:

\$\$v_C^c(S) = m(S_0, C^\*),\$\$

where \\m(S_0, C^\*)\\ is the cost of the minimum spanning tree of the
subproblem induced by the agents in \\S\\ and the source, evaluated over
the cycle-complete costs \\C^\*\\. The cycle-complete cost is defined as
the minimum between the direct cost \\c\_{ij}\\ and the smallest value
of the maximum arc cost over all cycles containing \\i\\ and \\j\\,
i.e., \\c\_{ij}^\* = \min(c\_{ij}, \min\_{f \in \mathcal{C}(i,j)}
\max\_{e \in f} c_e),\\ where \\\mathcal{C}(i,j)\\ denotes the set of
all cycles that include both \\i\\ and \\j\\.

For computational efficiency, \\C^\*\\ is obtained via the irreducible
matrices \\\bar{C}\\ of the subgraphs as:

\$\$c\_{ij}^\* = \max\_{k \in N \setminus \\i,j\\} {\bar{c}}\_{ij}^{N
\setminus \\k\\} \text{ for } i,j \in N,\$\$

and

\$\$c\_{0i}^\* = \max\_{k \in N \setminus \\i\\} {\bar{c}}\_{0i}^{N
\setminus \\k\\} \text{ for } i \in N,\$\$

where \\{\bar{C}}^{N \setminus \\k\\}\\ is the irreducible cost matrix
of the network excluding agent \\k\\. If no agents can be removed, the
irreducible cost of the full network is used.

Under `method = "exact"`, the function builds the full characteristic
function (all \\2^n - 1\\ values). Depending on the chosen solution
concept, it computes the allocation via
[`shapleyValue`](https://rdrr.io/pkg/CoopGame/man/shapleyValue.html)
(called *cycle-complete solution*, which relies on the coalition-based
formula) or
[`nucleolus`](https://rdrr.io/pkg/CoopGame/man/nucleolus.html).

For larger networks, `method = "montecarlo"` provides an approximation
of the Shapley value through random sampling using \\n\_{sim}\\
permutations. This approach implements lazy evaluation based on the
permutation formula of the Shapley value: instead of computing the MST
cost for all \\2^n - 1\\ possible subsets, it only evaluates the
specific subcoalitions generated by the random permutations.

This is highly useful when the number of agents is large, as it
significantly reduces memory usage and computation time while providing
a very accurate approximation of the *cycle-complete solution*. For
networks with few players, the exact method is preferred to avoid
potential overestimations.

For more details on these formulas, see
[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md).

## Note

The cycle-complete game is concave, which implies that its Shapley value
(the *cycle-complete solution*) belongs to its core. Moreover,
\\\mathrm{core}(N,v^c) \subseteq \mathrm{core}(N,v^p)\\.

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

Trudeau C (2012) A new stable and more responsible cost sharing solution
for mcst problems. Games Econom Behav 75(1):402–412

## See also

[`mcstGamePrivate`](https://luciasouto13.github.io/mcstprules/reference/mcstGamePrivate.md)
for the pessimistic game based on the original costs.

[`mcstGameIrred`](https://luciasouto13.github.io/mcstprules/reference/mcstGameIrred.md)
for the pessimistic game based on the irreducible cost matrix
\\\bar{C}\\.

[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md)
for an overview of the available rules and analysis tools in the
package.

## Examples

``` r
# Simple vector input
mcstGameCC(c(12, 15, 20, 4, 6, 8), sol = "shapley", draw = TRUE)


#>       S v^c(S)
#> 1   {1}     12
#> 2   {2}     15
#> 3   {3}     15
#> 4 {1,2}     16
#> 5 {1,3}     18
#> 6 {2,3}     23
#> 7     N     22
#> 
#> Solution: Shapley 
#>    1    2    3 
#> 4.33 8.33 9.33 

# Input with infinite costs (disconnected nodes)
C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
           Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
mcstGameCC(C_inf, sol = "shapley")
#>              S v^c(S)
#> 1          {1}      5
#> 2          {2}      8
#> 3          {3}      8
#> 4          {4}      8
#> 5          {5}      6
#> 6          {6}     11
#> 7        {1,2}     12
#> 8        {1,3}     13
#> 9        {1,4}     13
#> 10       {1,5}     11
#> 11       {1,6}     16
#> 12       {2,3}     16
#> 13       {2,4}     15
#> 14       {2,5}     14
#> 15       {2,6}     19
#> 16       {3,4}     13
#> 17       {3,5}     14
#> 18       {3,6}     19
#> 19       {4,5}     14
#> 20       {4,6}     17
#> 21       {5,6}     17
#> 22     {1,2,3}     20
#> 23     {1,2,4}     19
#> 24     {1,2,5}     18
#> 25     {1,2,6}     23
#> 26     {1,3,4}     18
#> 27     {1,3,5}     19
#> 28     {1,3,6}     24
#> 29     {1,4,5}     19
#> 30     {1,4,6}     22
#> 31     {1,5,6}     22
#> 32     {2,3,4}     20
#> 33     {2,3,5}     22
#> 34     {2,3,6}     27
#> 35     {2,4,5}     21
#> 36     {2,4,6}     24
#> 37     {2,5,6}     25
#> 38     {3,4,5}     19
#> 39     {3,4,6}     22
#> 40     {3,5,6}     25
#> 41     {4,5,6}     23
#> 42   {1,2,3,4}     24
#> 43   {1,2,3,5}     26
#> 44   {1,2,3,6}     31
#> 45   {1,2,4,5}     25
#> 46   {1,2,4,6}     28
#> 47   {1,2,5,6}     29
#> 48   {1,3,4,5}     24
#> 49   {1,3,4,6}     27
#> 50   {1,3,5,6}     30
#> 51   {1,4,5,6}     28
#> 52   {2,3,4,5}     26
#> 53   {2,3,4,6}     29
#> 54   {2,3,5,6}     33
#> 55   {2,4,5,6}     30
#> 56   {3,4,5,6}     28
#> 57 {1,2,3,4,5}     30
#> 58 {1,2,3,4,6}     33
#> 59 {1,2,3,5,6}     37
#> 60 {1,2,4,5,6}     34
#> 61 {1,3,4,5,6}     33
#> 62 {2,3,4,5,6}     35
#> 63           N     39
#> 
#> Solution: Shapley 
#>    1    2    3    4    5    6 
#>  4.5  7.0  6.5  5.0  6.0 10.0 

# Matrix input
C_mat <- matrix(c(0, 12, 15, 12,
                 12,  0,  4,  6,
                 15,  4,  0,  8,
                 12,  6,  8,  0), byrow = TRUE, ncol = 4)
mcstGameCC(C_mat, draw = TRUE, which.plot = "main", titles = FALSE)

#>       S v^c(S)
#> 1   {1}     12
#> 2   {2}     12
#> 3   {3}     12
#> 4 {1,2}     16
#> 5 {1,3}     18
#> 6 {2,3}     20
#> 7     N     22
#> 
#> Solution: Shapley 
#>    1    2    3 
#> 6.33 7.33 8.33 
```
