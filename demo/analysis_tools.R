# =====================================================================
#           mcstprules: Advanced Analysis & Stability Demo
# =====================================================================
# This demo explores the analytical tools available in the package.
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
costs <- data.frame(from = c(0, 0, 0, 1, 1, 2),
                    to   = c(1, 2, 3, 2, 3, 3),
                    cost = c(12, 15, 20, 4, 6, 8))
print(costs)
pause()

# ---------------------------------------------------------------------
#                      STEP 2: Building the Input
# ---------------------------------------------------------------------
# Some functions accept two types of input:

# Option A: derive the game from the cost matrix
priv <- mcstGamePrivate(costs)
print(priv)

# Option B: supply a characteristic function vector manually
vvals <- c(10, 12, 15, 20, 22, 25, 30)
vvals2 <- c(10, 10, 10, 15, 15, 15, 30)
pause()

# ---------------------------------------------------------------------
#                    STEP 3: Core Stability Check
# ---------------------------------------------------------------------
# mcstCoreCheck() verifies whether a given allocation belongs to the
# core.
check1 <- mcstCoreCheck(game = priv)
print(check1)

check2 <- mcstCoreCheck(allocation = c(15, 10, 5), v = vvals)
print(check2)
plot(check2)
pause()

# ---------------------------------------------------------------------
#                     STEP 4: Core Emptiness Check
# ---------------------------------------------------------------------
# mcstCorePoints() verifies whether the core is non-empty, returns its
# vertices and plot() renders the region.
cp1 <- mcstCorePoints(game = priv)
print(cp1)
plot(cp1)

cp2 <- mcstCorePoints(v = vvals2)
print(cp2)
pause()

# ---------------------------------------------------------------------
#                   STEP 5: Arc Sensitivity Analysis
# ---------------------------------------------------------------------
# How does an allocation rule respond to changes in a specific arc
# cost? mcstSensitivity() perturbs arc (i, j) by delta and tracks
# how each agent's payment evolves.
sens <- mcstSensitivity(costs, rule = "bird", arcs = c(0, 1))
plot(sens)

# =====================================================================
#            Demo completed. Explore more with ?mcstRules
# =====================================================================
