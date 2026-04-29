# =====================================================================
#                   mcstprules: Package Demonstration
# =====================================================================
# This demo shows how to compute and visualize cost allocation rules
# for Minimum Cost Spanning Tree Problems (MCSTP)
# =====================================================================

library(mcstprules)

# Helper function to pause execution for interactive use
pause <- function() {
  invisible(readline(prompt = "\n[Press Enter to continue...]"))
}

# ---------------------------------------------------------------------
#                     STEP 1: Defining the Problem
# ---------------------------------------------------------------------
# We define a cost matrix for a network with a source (node 0) and 3
# agents. Costs reflect the connection price between nodes.

costs <- matrix(c(0, 12, 15, 20,
                 12,  0,  4,  6,
                 15,  4,  0,  8,
                 20,  6,  8,  0), byrow = TRUE, ncol = 4)

print(costs)
pause()

# ---------------------------------------------------------------------
#                 STEP 2.1: Computing Algorithmic Rules
# ---------------------------------------------------------------------
# These rules are derived directly from the network structure and
# classical MST algorithms like Prim's and Kruskal's.

# 1. Bird's Rule: based on Prim's algorithm. Each agent pays the cost
#    of the arc that connects them to the source.
bird_res <- bird_rule(costs)

# 2. Folk Rule (ERO): based on Kruskal's algorithm. Agents divide the
#    cost of connecting their components equally.
folk_res <- folk_rule(costs)

# ---------------------------------------------------------------------
#               STEP 2.2: Computing Game-Theoretic Rules
# ---------------------------------------------------------------------
# These rules apply the Shapley value to cooperative games associated
# with the MCSTP.

# 1. Private Game (Kar's Rule): coalitions connect only through their
#    own members and the source (pessimistic approach).
kar_res <- private_game(costs)

# 2. Public Game: coalitions can connect through nodes outside the
#    coalition at no cost (optimistic approach).
public_res <- public_game(costs)

# ---------------------------------------------------------------------
#              STEP 3: Comparison of Allocations per Agent
# ---------------------------------------------------------------------
# Create a comparison table for all computed rules
summary_table <- rbind(
  "Bird"         = bird_res$bird,
  "Folk"         = folk_res$folk,
  "Private Game" = kar_res$shapley,
  "Public Game"  = public_res$shapley
)

print(round(summary_table, 2))
pause()

# ---------------------------------------------------------------------
#                   STEP 4: Visualizing the Networks
# ---------------------------------------------------------------------
# The package can plot the MST and how the total cost is shared among
# agents.

# Set layout
old_par <- par(mfrow = c(2,2), mar = c(2,2,4,2))

bird_rule(costs, draw = TRUE, which = 2)
folk_rule(costs, draw = TRUE)
private_game(costs, draw = TRUE, which = 2)
public_game(costs, draw = TRUE, which = 2)

par(old_par)
pause()

# ---------------------------------------------------------------------
#                         STEP 5: Handling Ties
# ---------------------------------------------------------------------
# When multiple MCSTs exist, rules based on Prim's algorithm can vary.
# The package automatically computes the symmetric average over all
# permutations (Extended Bird Rule).

# Symmetric case: agents 1 and 2 are at the same distance from the
# source
tie_costs <- c(10, 10, 2) # Format: (0,1), (0,2), (1,2)
bird_tie <- bird_rule(tie_costs)

print(bird_tie$perms)
print(bird_tie$e_bird)

# Note: Kruskal-based rules like the folk rule are unaffected by ties.

# =====================================================================
#            Demo completed. Explore more with ?alloc_rules
# =====================================================================
