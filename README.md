# pyDualBB

[![CodeFactor](https://www.codefactor.io/repository/github/jhelgert/pydualbb/badge)](https://www.codefactor.io/repository/github/jhelgert/pydualbb)

This python module provides python bindings for a so-called *dual Branch-and-Bound*
(dualBB) algorithm in order to solve convex mixed-integer quadratically
constrained quadratic programs of the form:

```
min  f(x) = 0.5 * x' * Q0 * x + c0' * x + r0
 x
s.t. g_i(x) = 0.5 * x' * Qi * x + ci' * x + r0 <= 0 for all i = 1, .., p
      h1(x) = A1*x - b1 <= 0
      h2(x) = A2*x - b2 == 0
      x_i integer for all i in I
```

for a small number `p`, i.e. a small number of quadratic constraints.

Here,

- all the matrices `Q0`, ... `Qp` are assumed to be symmetric positive definite,
i.e. `f` and `g_i` are strictly convex.
- the matrices `A1` and `A2` have dimensions `(m1, n)` and `(m2,n)` with `m2 <= n` and `A2` has full rank.

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
