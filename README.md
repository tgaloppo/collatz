# Description
Companion code to: 

"Negative Drift and State Instability in a Bitwise System Equivalent to the Collatz Conjecture"

Preprint is available here: https://doi.org/10.5281/zenodo.17407872

and

"Deterministic Limits and Ergodic Properties of a Bitwise Syracuse Map"

Preprint is available here: https://doi.org/10.5281/zenodo.17544884

This repository contains code and scripts to reproduce results in the above research papers regarding the Stochastic Stability of the Collatz Conjecture. It implements a bitwise map simulation and statistical analyses of bit-span dynamics, and a Ulam Method estimatation of the transfer operator's spectral gap.

# Requirements

This code requires:

* The GNU MP BigNum library (GMP), available here: https://gmplib.org/
* The Eigen3 library for linear algebra, available here: https://eigen.tuxfamily.org
* The Spectra library for large scale eigenvalue problems, available here: https://spectralib.org/

If you are using a Ubuntu based system you should be able to install all dependencies with apt:

$ sudo apt install libeigen3-dev libspectra-dev libgmp-dev

# Building

The supplied Makefile should build everything for you:

$ make

# Running

## Main program

The primary executable that computes the bit-span metrics can accept the number of random starting values and the maximum starting value size, in bits:

**collatz** [num values] [num bits]

The defaults are "num values" = 1000000 and "num bits" = 128.

## Recreating the plots

If you want to recreate the plots from the papers, or to generate the data for your own exploration, use

$ make generate-plots

or

$ make generate-data

respectively. Expect it to take about 20 minutes. Note that plot generation requires R with ggplot2 and dplyr packages.

Generated plots include:
* **benford-dist.eps** - Plots the distribution of the fractional position of iterates (Paper 2, Figure 2)
  * The data for this plot is **fractional-position.csv**
* **autocorrelation.eps** - Plot the autocorrelation of the fractional positions over sequence lags (Paper 2, Figure 3)
  * The data for this plot is **autocorrelation.csv**
* **spectral-gap.eps** - Log-log plot of the spectral gap vs input bit size (Supplement, Figure 1)
  * The data for this plot is **ulam-gap-out.csv**
* **delta-L-freq.eps** - Semi-log plot of LSB shift distribution (Supplement, Figure 2)
  * The data for this plot is **delta-L-freq.csv**
