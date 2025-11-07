# Description
Companion code to: 

"Negative Drift and State Instability in a Bitwise System Equivalent to the Collatz Conjecture"
Preprint is available here: https://zenodo.org/records/17407873

and

"Deterministic Limits and Ergodic Properties of a Bitwise Syracuse Map"
Preprint is available here: https://zenodo.org/records/17552420

This code reproduces all numerical results presented in the papers.

# Requirements

This code requires:

* The GNU MP BigNum library (GMP), available here: https://gmplib.org/
* The Eigen3 library for linear algebra, available here: https://eigen.tuxfamily.org

# Building

## Typical

$ g++ -o collatz collatz.cpp -lgmp -O3

## If your Eigen headers are not in standard place:

$ g++ -o collatz collatz.cpp -lgmp -O3 -I /path/to/eigen3

# Running

The executable can accept the number of random starting values and the maximum starting value size, in bits:

**collatz** [num values] [num bits]

The defaults are "num values" = 1000000 and "num bits" = 128.
