# pyDualBB

This python module provides python bindings for a so-called *dual Branch-and-Bound*
(dualBB) algorithm in order to solve convex mixed-integer quadratically
constrained quadratic programs of the form:

```math
\begin{align*}
\min_{x} f(x) &= \frac12 x^\top Q_0 x + c_0^\top x + r_0 \\
\text{s.t.} \quad g_k(x) &= \frac12 x^\top Q_k x + c_k^\top x + r_k \leq 0 \quad \forall k = 1,\ldots, p \\
      h_1(x) &= A_1 x - b_1 \leq 0 \\
      h_2(x) &= A_2 x - b_2 = 0 \\
        x_i &\in \mathbb{Z} \quad \forall i \in I
\end{align*}
```

for a small number $p$ of quadratic constraints. Note that,
- all the matrices $Q_0$, ... $Q_p$ are assumed to be symmetric positive definite, which implies that the functions $f, g_1, \ldots, g_p$ are strictly convex.
- $A_1 \in \mathbb{R}^{m_1 \times n}$, $A_2 \in \mathbb{R}^{m_2 \times n}$ and $m_2 \leq n$.
- the matrix $A_2$ has full rank.

The algorithm is a multi-threaded and high-performance version of the original
implementation. See [here](https://link.springer.com/chapter/10.1007%2F978-3-030-48439-2_15) for a detailed description of the algorithm.

It is written in C++ and builds on CPLEX'
*generic callbacks*. The packages requires an installation of Cmake (>= 3.16),
a decent C++ compiler supporting C++17 and an installation of Cplex 20.10.


## Install

Clone this repo and run

```
python3 setup.py install
```
inside the repo folder.

## Example

``` python
import numpy as np
from pyDualBB import dualBB

# Set the cplex options for dualBB
# See the Cplex documentation for the specific values
cplex_options_dualbb = {
    "mip_strategy_nodeselect": 2,       # best estimate search
    "mip_strategy_search": 0,           # let cplex choose (default)
    "mip_strategy_heuristiceffort": 1,  # default
    "timelimit": 1800                   # 30 min
}

# Path to the problem as .lp / .mps or .mps.gz file
problem_file = "problem.mps.gz"
# Path to the .json file containing all matrices and vectors
json_file = "problem.json"

# Fixed lagrangian multipliers
alphas = np.array([0.01, 0.025, 0.05, 0.1, 0.2, 0.95])

# set the number of threads used for dualbb (all other available threads will be used by Cplex)
dualbb_threads = 4

status, total_runtime, incumbents, incRuntimes = dualBB(problem_file, json_file, alphas, options=cplex_options_dualbb, dual_branching_threads=dualbb_threads)

# Always cut-off the branch with higher dual bound (only for demonstration purposes!)
status, total_runtime, incumbents, incRuntimes = dualBB(problem_file, json_file, alphas, options=cplex_options_dualbb, only_single_branch=True, dual_branching_threads=dualbb_threads)
```

## To do:

- Provide a matrix interface for dualBB instead of using .mps and .json files
- Parse all model matrices directly from the .mps/.lp/.mps.gz file
