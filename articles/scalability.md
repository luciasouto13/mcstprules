# Demo for 'mcstprules::scalability'

This demo runs the code from
`demo("scalability", package = "mcstprules")` in R.

``` r

# =====================================================================
#              mcstprules: Scalability & Monte Carlo Demo
# =====================================================================
# This demo illustrates how mcstprules handles computationally
# intensive scenarios: large networks with multiple optimal trees and
# large cooperative games using Monte Carlo simulation.
# =====================================================================

library(mcstprules)

# Helper function to pause execution for interactive use
pause <- function() {
  invisible(readline(prompt = "\n[Press Enter to continue...]"))
}

# ---------------------------------------------------------------------
#              STEP 1: Defining a Large Problem with Ties
# ---------------------------------------------------------------------
# We generate a random symmetric cost matrix for 20 agents (+ source).
# An exact calculation here would require 20! permutations.
set.seed(123)
N0 <- 21 # 1 source (node 0) + 20 agents (nodes 1-20)
costs <- matrix(0, nrow = N0, ncol = N0)

# Cheap agent-to-agent connections
costs[lower.tri(costs)] <- sample(1:3, (N0 * (N0 - 1)) / 2,
                                  replace = TRUE)
costs <- costs + t(costs)

# Expensive source-to-agent connections
costs[1, 2:N0] <- sample(10:12, N0 - 1, replace = TRUE)
costs[2:N0, 1] <- sample(10:12, N0 - 1, replace = TRUE)
diag(costs) <- 0

print(head(costs[, 1:6]))
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0   11   10   11   10   12
#> [2,]   11    0    1    1    3    2
#> [3,]   12    1    0    1    2    3
#> [4,]   10    1    1    0    3    1
#> [5,]   11    3    2    3    0    2
#> [6,]   11    2    3    1    2    0
pause()
#> 
#> [Press Enter to continue...]

# ---------------------------------------------------------------------
#            STEP 2: Approximating Rules with Multiple MCSTs
# ---------------------------------------------------------------------
# When multiple minimal trees exist, rules based on tree-building
# procedures requires computing the symmetric average over
# permutations. At this scale (20!), exact enumeration is impossible.
# We use Monte Carlo sampling with Bird's rule as an example.
bs <- Sys.time()
bird <- mcstBird(costs, method = "montecarlo", nsim = 2000)
be <- Sys.time()

time_bird <- round(as.numeric(difftime(be, bs,
                                       units = "secs")), 4)
print(time_bird)
#> [1] 3.4733
pause()
#> 
#> [Press Enter to continue...]

# ---------------------------------------------------------------------
#          STEP 3: Approximating Game-Theoretic Rules at Scale
# ---------------------------------------------------------------------
# Similarly, computing the exact Shapley value for 20 agents is
# computationally impossible. The Monte Carlo method samples
# permutations to solve the game efficiently, exemplified here with
# the Private Game (Kar's allocation).
ks <- Sys.time()
kar <- mcstGamePrivate(costs, method = "montecarlo", nsim = 2000)
ke <- Sys.time()

time_kar <- round(as.numeric(difftime(ke, ks,
                                     units = "secs")), 4)
print(time_kar)
#> [1] 6.6336
pause()
#> 
#> [Press Enter to continue...]

# ---------------------------------------------------------------------
#                     STEP 4: Viewing the Solutions
# ---------------------------------------------------------------------
results <- data.frame(
  Agent    = 1:20,
  Bird_MC  = round(bird$e_bird, 2),
  Kar_MC   = round(kar$allocation, 2)
)
print(head(results))
#>   Agent Bird_MC Kar_MC
#> 1     1    1.00   1.48
#> 2     2    2.10   1.23
#> 3     3    1.00   1.49
#> 4     4    2.04   1.33
#> 5     5    1.00   1.70
#> 6     6    1.00   1.75

# =====================================================================
#             Demo completed. Explore more with ?mcstRules
# =====================================================================
```

------------------------------------------------------------------------

You can also [download the raw R script
here](https://github.com/luciasouto13/mcstprules/blob/master/demo/scalability.R)
to run it locally on your computer.
