# Cost Allocations based on the Private Game for MCSTP

This function computes the cost allocation for the private or
pessimistic game (Bird, 1976) associated with a minimum cost spanning
tree problem \\(N_0, C)\\, using a specified solution concept. The
private game, denoted as \\(N, v_C^p)\\ \\(\\or simply \\(N,
v^p)\\\\)\\, evaluates the cost of each coalition \\S \subseteq N\\
under the assumption that agents outside the coalition are not available
to provide their nodes as intermediaries.

## Usage

``` r
mcstGamePrivate(
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

- `v_p`: a numeric vector containing the characteristic function values
  \\v^p(S)\\ for all coalitions. If `method = "montecarlo"`, returns a
  string indicating that lazy evaluation was used.

- `contributions`: a data frame with the marginal contributions for the
  sampled permutations. If `method = "exact"`, returns a string
  indicating that this table is not used.

- `allocation`: a numeric vector containing the computed cost allocation
  for each agent.

- `total`: the total cost of the MCST, \\v^p(N)\\.

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

The private game \\(N, v_C^p)\\ is a cooperative game with transferable
utility (TU game) defined for each coalition \\S \subseteq N\\ as:

\$\$v_C^p(S) = m(S_0, C),\$\$

where \\m(S_0, C)\\ is the cost of the minimum spanning tree of the
subproblem induced by the agents in \\S\\ and the source. This approach
is termed "private" or "pessimistic" because the nodes in \\N \setminus
S\\ belong to these agents, and their participation is needed in order
to use their nodes.

When the Shapley value (`sol = "shapley"`) is applied to this game, the
resulting allocation is known as Kar's rule (Kar, 2002). For a function
that particularizes this specific case, see
[`mcstKar`](https://luciasouto13.github.io/mcstprules/reference/mcstKar.md).

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
a very accurate approximation of Kar's rule. For networks with few
players, the exact method is preferred to avoid potential
overestimations.

For more details on these formulas, see
[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md).

## Note

The private game is a strong reference for stability. Classical rules
such as Bird, Dutta-Kar, or the folk rule are known to belong to the
core of \\v^p\\ for any MCST problem. As a result, these allocations are
coalitionally stable under the pessimistic cost perspective.

Although the Shapley value and the nucleolus are well-established
solution concepts, it has been proven that computing both for this game
is NP-hard.

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

Bird C (1976) On cost allocation for a spanning tree: a game theoretic
approach. Networks 6(4):335–350

Kar A (2002) Axiomatization of the Shapley value on minimum cost
spanning tree games. Games Econom Behav 38(2):265–277

## See also

[`mcstBird`](https://luciasouto13.github.io/mcstprules/reference/mcstBird.md),
[`mcstDuttaKar`](https://luciasouto13.github.io/mcstprules/reference/mcstDuttaKar.md),
[`mcstFolk`](https://luciasouto13.github.io/mcstprules/reference/mcstFolk.md)
for classical rules belonging to the core of \\v^p\\.

[`mcstKar`](https://luciasouto13.github.io/mcstprules/reference/mcstKar.md)
for the specific cost allocation rule based on the Shapley value of this
game.

[`mcstGameIrred`](https://luciasouto13.github.io/mcstprules/reference/mcstGameIrred.md),
[`mcstGameCC`](https://luciasouto13.github.io/mcstprules/reference/mcstGameCC.md)
for the pessimistic games based on the irreducible cost matrix
\\\bar{C}\\ and the cycle-complete cost matrix \\C^\*\\, respectively.

[`mcstGameOpt`](https://luciasouto13.github.io/mcstprules/reference/mcstGameOpt.md),
[`mcstGamePublic`](https://luciasouto13.github.io/mcstprules/reference/mcstGamePublic.md)
for optimistic alternatives to the private approach.

[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md)
for an overview of the available rules and analysis tools in the
package.

## Examples

``` r
# Simple vector input
mcstGamePrivate(c(12, 15, 20, 4, 6, 8), sol = "shapley", draw = TRUE)


#>       S v^p(S)
#> 1   {1}     12
#> 2   {2}     15
#> 3   {3}     20
#> 4 {1,2}     16
#> 5 {1,3}     18
#> 6 {2,3}     23
#> 7     N     22
#> 
#> Solution: Shapley 
#>    1    2    3 
#>  3.5  7.5 11.0 

# Input with infinite costs (disconnected nodes)
C_inf <- c(5, 9, Inf, 15, 6, Inf, 7, Inf, Inf, Inf,
           Inf, 8, 7, Inf, Inf, 5, Inf, Inf, 8, 9, 11)
mcstGamePrivate(C_inf, sol = "nucleolus")
#>              S v^p(S)
#> 1          {1}      5
#> 2          {2}      9
#> 3          {3}      0
#> 4          {4}     15
#> 5          {5}      6
#> 6          {6}      0
#> 7        {1,2}     12
#> 8        {1,3}      5
#> 9        {1,4}     20
#> 10       {1,5}     11
#> 11       {1,6}      5
#> 12       {2,3}     17
#> 13       {2,4}     16
#> 14       {2,5}     15
#> 15       {2,6}      9
#> 16       {3,4}     20
#> 17       {3,5}      6
#> 18       {3,6}      0
#> 19       {4,5}     14
#> 20       {4,6}     24
#> 21       {5,6}     17
#> 22     {1,2,3}     20
#> 23     {1,2,4}     19
#> 24     {1,2,5}     18
#> 25     {1,2,6}     12
#> 26     {1,3,4}     25
#> 27     {1,3,5}     11
#> 28     {1,3,6}      5
#> 29     {1,4,5}     19
#> 30     {1,4,6}     29
#> 31     {1,5,6}     22
#> 32     {2,3,4}     21
#> 33     {2,3,5}     23
#> 34     {2,3,6}     17
#> 35     {2,4,5}     21
#> 36     {2,4,6}     25
#> 37     {2,5,6}     26
#> 38     {3,4,5}     19
#> 39     {3,4,6}     29
#> 40     {3,5,6}     17
#> 41     {4,5,6}     23
#> 42   {1,2,3,4}     24
#> 43   {1,2,3,5}     26
#> 44   {1,2,3,6}     20
#> 45   {1,2,4,5}     25
#> 46   {1,2,4,6}     28
#> 47   {1,2,5,6}     29
#> 48   {1,3,4,5}     24
#> 49   {1,3,4,6}     34
#> 50   {1,3,5,6}     22
#> 51   {1,4,5,6}     28
#> 52   {2,3,4,5}     26
#> 53   {2,3,4,6}     30
#> 54   {2,3,5,6}     34
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
#> Solution: Nucleolus 
#>  1  2  3  4  5  6 
#>  5  9 -2 15  6  6 

# Matrix input
C_mat <- matrix(c(0, 12, 15, 12,
                 12,  0,  4,  6,
                 15,  4,  0,  8,
                 12,  6,  8,  0), byrow = TRUE, ncol = 4)
mcstGamePrivate(C_mat, draw = TRUE, which.plot = "main", titles = FALSE)

#>       S v^p(S)
#> 1   {1}     12
#> 2   {2}     15
#> 3   {3}     12
#> 4 {1,2}     16
#> 5 {1,3}     18
#> 6 {2,3}     20
#> 7     N     22
#> 
#> Solution: Shapley 
#>    1    2    3 
#> 5.83 8.33 7.83 
```
