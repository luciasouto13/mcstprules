# Stability Range for MCST Arcs

This function calculates the stability range for the arcs of a minimum
cost spanning tree (MCST). It determines how much the cost of an
individual arc can increase or decrease without changing the topology of
the optimal tree.

## Usage

``` r
mcstStabilityRange(C)
```

## Arguments

- C:

  cost matrix between nodes. Accepts multiple formats (see "Supported
  formats for `C`" in
  [`mcstBird`](https://luciasouto13.github.io/mcstprules/reference/mcstBird.md)).

## Value

A data frame containing:

- `from`: the starting node of the arc.

- `to`: the ending node of the arc.

- `in_mst`: logical; if `TRUE`, the arc belongs to the optimal tree.

- `cost`: the original cost of the arc \\c\_{ij}\\.

- `range`: a formatted string showing the allowed cost variation
  \\\small \[\Delta LB, \Delta UB\]\\.

- `LB`: the absolute lower bound for the arc's cost.

- `UB`: the absolute upper bound for the arc's cost.

## Details

The stability range of an arc \\(i, j)\\ indicates the interval \\\[LB,
UB\]\\ within which its original cost \\c\_{ij}\\ can vary while
preserving the current minimum cost spanning tree as optimal.

At the finite boundary points of the stability interval, the minimum
cost spanning tree becomes non-unique due to cost ties. In these
scenarios, component or path-based rules (e.g., Folk, Kar) remain
stable, whereas permutation-based rules such as Bird's rule or
Dutta-Kar's rule average allocations over all alternative optimal trees,
which can induce non-proportional jumps in the final cost allocation
vector.

## See also

[`mcstBird`](https://luciasouto13.github.io/mcstprules/reference/mcstBird.md),
[`mcstDuttaKar`](https://luciasouto13.github.io/mcstprules/reference/mcstDuttaKar.md),
[`mcstFolk`](https://luciasouto13.github.io/mcstprules/reference/mcstFolk.md),
[`mcstKar`](https://luciasouto13.github.io/mcstprules/reference/mcstKar.md)
for allocation rules.

[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md)
for an overview of the available rules and analysis tools in the
package.

## Examples

``` r
# Matrix input
C_mat <- matrix(c(0, 10, 15, 20,
                 10,  0, 25, 12,
                 15, 25,  0,  8,
                 20, 12,  8,  0), byrow = TRUE, ncol = 4)

# Calculate stability range
mcstStabilityRange(C_mat)
#>  from to in_mst cost      range LB  UB
#>     0  1   TRUE   10   [-10, 5]  0  15
#>     0  2  FALSE   15  [-3, Inf] 12 Inf
#>     0  3  FALSE   20  [-8, Inf] 12 Inf
#>     1  2  FALSE   25 [-13, Inf] 12 Inf
#>     1  3   TRUE   12   [-12, 3]  0  15
#>     2  3   TRUE    8    [-8, 7]  0  15
```
