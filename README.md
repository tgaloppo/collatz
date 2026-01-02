# Description
Companion code to: 

[Negative Drift and State Instability in a Bitwise System Equivalent to the Collatz Conjecture][1]

and

[Deterministic Limits and Ergodic Properties of a Bitwise Syracuse Map][2]

This repository contains code and scripts to reproduce results in the above research papers regarding the Stochastic Stability of the Collatz Conjecture.
It implements an arbitrary-precision bitwise map simulation capable of analyzing trajectories exceeding 10^7 bits and statistical analyses of bit-span dynamics, and a Ulam Method estimation of the transfer operator's spectral gap.

## Mathematical Significance: The Zero-Noise Limit

This framework addresses the Collatz Conjecture by reducing it to a specific problem of **stochastic stability** in the zero-noise limit.

### The Reduction
The code implements a bitwise dynamical system $U(n)$ which acts as a perturbed irrational rotation on the 2-adic integers. We define the dynamics as:
$$T(x) = T_{asymp}(x) + \epsilon_j(x)$$
Where $T_{asymp}$ is a pure rotation (ergodic but non-mixing) and $\epsilon_j$ is a state-dependent perturbation that decays as the bit-span $S \to \infty$.

### The Open Question
The central contribution of this work is the formal reduction of the conjecture to a single open question of ergodicity: **Does the invariant measure of the system remain absolutely continuous in the limit as $\|\epsilon\| \to 0$?**

* **ACIM Stability:** We prove that if the invariant measure possesses *any* continuous probability density function (not necessarily Benford), the low-order bits of the orbit must be uniformly distributed.
* **Negative Drift:** A uniform distribution of low-order bits necessitates a global negative drift in the bit-span metric ($E[\Delta S] \approx -0.41$), rendering divergence impossible.
* **Singular Collapse:** Divergence is therefore only possible if the invariant measure collapses onto a singular support as the perturbation vanishes.

### Empirical Validation
This repository provides the tools to test this stability. The `collatz` executable can simulate trajectories exceeding **$10^7$ bits**, demonstrating that the system's mixing properties and drift metrics remain scale-invariant even as the perturbation magnitude approaches $2^{-10,000,000}$. This suggests the system is stochastically stable and the absolutely continuous measure persists.

# Requirements

This code requires:

* The GNU MP BigNum library (GMP), available here: https://gmplib.org/
* The Eigen3 library for linear algebra, available here: https://eigen.tuxfamily.org
* The Spectra library for large scale eigenvalue problems, available here: https://spectralib.org/

If you are using a Ubuntu based system you should be able to install all dependencies with apt:

```sudo apt install libeigen3-dev libspectra-dev libgmp-dev```

# Building

The supplied Makefile should build everything for you:

```make```

# Running

## Main program

The primary executable that computes the bit-span metrics can accept the number of random starting values and the maximum starting value size, in bits:

```./collatz [num values] [num bits]```

The defaults are "num values" = 1000000 and "num bits" = 128. You can run singular trajectories with massive starting values such as:

```./collatz 1 1000000```

## Recreating the plots

If you want to recreate the plots from the papers, or to generate the data for your own exploration, use

```make generate-plots```

or

```make generate-data```

respectively. Expect it to take about 20 minutes. Note that plot generation requires R with ggplot2 and dplyr packages. Installing these
on Ubuntu based systems can be accomplished with:

```sudo apt install r-base r-cran-ggplot2 r-cran-dplyr```

Generated plots include:
* **benford-dist.eps** - Plots the distribution of the fractional position of iterates (Paper 2, Figure 2)
  * The data for this plot is **fractional-position.csv**
* **autocorrelation.eps** - Plot the autocorrelation of the fractional positions over sequence lags (Paper 2, Figure 3)
  * The data for this plot is **autocorrelation.csv**
* **spectral-gap.eps** - Log-log plot of the spectral gap vs input bit size (Supplement, Figure 1)
  * The data for this plot is **ulam-gap-out.csv**
* **delta-L-freq.eps** - Semi-log plot of LSB shift distribution (Supplement, Figure 2)
  * The data for this plot is **delta-L-freq.csv**

[1]: https://doi.org/10.5281/zenodo.17407872
[2]: https://doi.org/10.5281/zenodo.17544884
