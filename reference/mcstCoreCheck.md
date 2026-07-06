# Core Stability Check for Cost Allocations

This function evaluates whether a given cost allocation belongs to the
*core* of a cooperative game. It checks for both efficiency and
coalitional rationality, identifying any blocking coalitions if the
allocation is unstable.

## Usage

``` r
mcstCoreCheck(game = NULL, allocation = NULL, v = NULL, tol = 1e-07)

# S3 method for class 'mcstp_core'
plot(x, col = "firebrick", titles = TRUE, ...)
```

## Arguments

- game:

  list; an object containing a cooperative game computed with
  `method = "exact"` (e.g., from
  [`mcstGamePrivate`](https://luciasouto13.github.io/mcstprules/reference/mcstGamePrivate.md)).
  Default is `NULL`.

- allocation:

  numeric vector; a proposed cost allocation for each agent. Required if
  `game` is `NULL`.

- v:

  numeric vector; the characteristic function values \\v(S)\\ for all
  \\2^n - 1\\ coalitions, ordered by cardinality. Required if `game` is
  `NULL`.

- tol:

  numeric; tolerance for floating-point comparisons to avoid precision
  issues. Default is `1e-7`.

- x:

  the output from `mcstCoreCheck`.

- col:

  character string specifying the color for the excess bars. Default is
  `"firebrick"`.

- titles:

  logical; if `TRUE` (default), adds a main title and a y-axis label to
  the plot.

- ...:

  additional graphical parameters passed to
  [`barplot`](https://rdrr.io/r/graphics/barplot.html).

## Value

A list containing:

- `in_core`: logical; if `TRUE`, the allocation is stable (belongs to
  the *core*).

- `is_efficient`: logical; if `TRUE`, the sum of the allocation equals
  the total cost.

- `is_rational`: logical; if `TRUE`, no blocking coalitions exist.

- `blocking_coals`: a data frame detailing the blocking coalitions
  \\S\\, their allocated cost \\x(S)\\, their value \\v(S)\\, and the
  excess. Empty if the allocation is rational.

- `n`: the number of agents.

- `sum_x`: the sum of the allocated costs.

- `total`: the total cost of the game \\v(N)\\.

## Details

An allocation \\x \in \mathbb{R}^N\\ belongs to the *core* of a
cooperative game \\(N, v)\\ if it satisfies two fundamental conditions:

1.  Efficiency: the total cost is completely allocated among the agents,
    i.e., \$\$\sum\_{i \in N} x_i = v(N).\$\$

2.  Coalitional rationality: no coalition pays more than its stand-alone
    cost, i.e., \$\$\sum\_{i \in S} x_i \le v(S) \quad \forall S \subset
    N.\$\$

If the coalitional rationality condition is violated for any coalition
\\S\\, \\S\\ is considered a *blocking coalition*, as its members would
have an incentive to leave the grand coalition and form their own
sub-network. The excess \\x(S) - v(S)\\ quantifies the magnitude of this
violation.

The function can be executed in two ways: by passing a complete `game`
object (which must have been computed with `method = "exact"`), or by
manually providing the proposed `allocation` vector and the
characteristic function `v`.

The `plot` method visualizes the excess of the blocking coalitions (if
any) using a bar chart.

## References

Bergantiños G, Vidal-Puga J (2021) A review of cooperative rules and
their associated algorithms for minimum-cost spanning tree problems.
SERIEs, 12:73-100

## See also

[`mcstGamePrivate`](https://luciasouto13.github.io/mcstprules/reference/mcstGamePrivate.md),
[`mcstGameIrred`](https://luciasouto13.github.io/mcstprules/reference/mcstGameIrred.md)
for computing cooperative games that can be passed to this function.

[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md)
for an overview of the available rules and analysis tools in the
package.

## Examples

``` r
# Example 1: Using manual vectors
vvals <- c(10, 12, 15, 20, 22, 25, 30) # 3 agents: v(1), v(2), v(3), v(12), ...
# Stable check
alloc1 <- c(8, 10, 12)
check1 <- mcstCoreCheck(allocation = alloc1, v = vvals); check1
#> -------------------------
#>  Core Stability Analysis
#> -------------------------
#> Players (n) : 3 
#> Efficiency  : Passed (Sum x: 30.00 | v(N): 30.00)
#> Rationality : Passed 
#> 
#> Result: the allocation IS IN THE CORE (Stable)
# Unstable check
alloc2 <- c(15, 10, 5)
check2 <- mcstCoreCheck(allocation = alloc2, v = vvals); check2
#> -------------------------
#>  Core Stability Analysis
#> -------------------------
#> Players (n) : 3 
#> Efficiency  : Passed (Sum x: 30.00 | v(N): 30.00)
#> Rationality : Failed 
#> 
#> Result: the allocation IS NOT IN THE CORE (Unstable)
#> 
#> Blocking Coalitions (x(S) > v(S)):
#>      S x(S) v(S) excess
#>    {1}   15   10      5
#>  {1,2}   25   20      5
plot(check2)

# \donttest{
# Example 2: Using a game object
C <- matrix(c(0, 12, 15, 12, 0, 4, 15, 4, 0), byrow = TRUE, ncol = 3)
# Check if the Shapley value of the private game is in the core
mcstCoreCheck(game = mcstGamePrivate(C))
#> -------------------------
#>  Core Stability Analysis
#> -------------------------
#> Players (n) : 2 
#> Efficiency  : Passed (Sum x: 16.00 | v(N): 16.00)
#> Rationality : Passed 
#> 
#> Result: the allocation IS IN THE CORE (Stable)
# }
```
