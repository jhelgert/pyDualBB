# pyDualBB

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

##