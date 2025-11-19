// This code is a companion to the paper "Deterministic Limits and Ergodic
// Properties of a Bitwise Syracuse Map". This program uses Ulam's method
// to estimate the spectral gap of the transfer operator of the map U(a)
// as defined in the paper. The spectral gap is estimated increasing
// bit sizes for comparison with the vanishing noise term (epsilon).
//
// This code requires:
//   1. The GNU MP BigNum library (GMP), available here: https://gmplib.org/
//   2. The Eigen3 library for linear algebra, available here: https://eigen.tuxfamily.org
//   3. The Spectra library for large scale eigenvalue problems, available here: https://spectralib.org/
//
// To build:
// g++ -o ulam-gap ulam-gap.cpp -lgmp -O3
//
// You may need to specify the location of your Eigen3 headers:
// g++ -o ulam-gap ulam-gap.cpp -lgmp -O3 -I /path/to/eigen3
//
// Auth: Travis Galoppo (tgaloppo@gmail.com)
// Date: 2025-11-19
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <gmp.h>
#include <eigen3/Eigen/Sparse>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

constexpr size_t NUM_BINS = 10000;

void ulam(size_t num_values, size_t num_bits, Eigen::SparseMatrix<double>* P) {
  mpz_t n, nu, two, tmp;
  gmp_randstate_t rnd_state;
  long exp_temp;
  
  // initialize
  gmp_randinit_default(rnd_state);
  gmp_randseed_ui(rnd_state, time(NULL));
  
  mpz_inits(n, nu, two, tmp, NULL);
  
  mpz_set_ui(two, 2);

  size_t L0, L1;
  size_t H0, H1;

  // main loop
  for (size_t j = 0; j < num_values; j++) {
    // generate a random starting value
    mpz_urandomb(n, rnd_state, num_bits);
  
    // find LSB, MSB and initial span
    L0 = mpz_scan1(n, 0);
    H0 = mpz_sizeinbase(n, 2) - 1;

    // determine initial bin
    double alpha = mpz_get_d_2exp(&exp_temp, n);
    size_t prev_bin = (size_t)((alpha - 0.5) * 2.0 * NUM_BINS);
    // Safety clamp (floating point errors can sometimes hit 1.0 exactly)
    if (prev_bin >= NUM_BINS) prev_bin = NUM_BINS - 1;  

    // loop until we reach a power of 2
    while ( H0 > L0 ) {
      // 2-adic value of n
      mpz_pow_ui(nu, two, L0);

      // apply U(n)
      mpz_mul_ui(n, n, 3);
      mpz_add(n, n, nu);
    
      // get fractional position of new value
      alpha = mpz_get_d_2exp(&exp_temp, n);

      // Map alpha [0.5, 1.0) -> Bin Index [0, NUM_BINS-1]
      size_t curr_bin = (size_t)((alpha - 0.5) * 2.0 * NUM_BINS);
      // Safety clamp (floating point errors can sometimes hit 1.0 exactly)
      if (curr_bin >= NUM_BINS) curr_bin = NUM_BINS - 1;

      // Record transition from previous step
      P->coeffRef(prev_bin, curr_bin) += 1.0;

      prev_bin = curr_bin;

      // find LSB, MSB, and span of resulting number
      L1 = mpz_scan1(n, L0); // L1 must be >L0, so this safe
      H1 = mpz_sizeinbase(n,2) - 1;
   
      // Update starting conditions for next iteration;
      // for efficiency, pull LSB back to position 0
      mpz_tdiv_q_2exp(n, n, L1);
      L0 = 0;
      H0 = H1 - L1;
    }
  }
}

int main(int argc, char** argv) {
  size_t ncpu = std::thread::hardware_concurrency();

  // default number of values and bit size
  size_t num_values = 1000000;
  size_t num_bits   = 100;

  while (num_bits < 10000000) {
    std::cerr << "Computing Matrix @ " << num_bits << " bits" << std::endl;

    std::vector<Eigen::SparseMatrix<double>> P(ncpu);
    for (size_t j = 0; j < ncpu; j++) {
        P[j].resize(NUM_BINS, NUM_BINS);
    }

    size_t n_per_thread = num_values / ncpu;

    std::vector<std::thread> threads(ncpu-1);
    for (size_t j = 1; j < ncpu; j++) {
        threads[j-1] = std::thread(ulam, n_per_thread, num_bits, &P[j]);
    }

    // main thread potentially gets to do a few extra...
    size_t n_final_thread = num_values - (ncpu - 1) * n_per_thread;
    ulam(n_final_thread, num_bits, &P[0]);

    // sum all partial matrices into P[0]
    for (size_t j = 0; j < threads.size(); j++ ) {
        threads[j].join();
        P[0] += P[j+1];
    }

    // Normalize rows of matrix
    for (size_t j = 0; j < NUM_BINS; j++) {
      double sum = P[0].row(j).sum();
      if (sum > 0.0) 
        P[0].row(j) /= sum;
    }

    std::cerr << "Computing Spectral Gap @ " << num_bits << " bits" << std::endl;

    // Construct matrix operation object using the wrapper class SparseGenMatProd
    Spectra::SparseGenMatProd<double> op(P[0]);

    // Construct eigen solver object, requesting the largest three eigenvalues
    Spectra::GenEigsSolver<Spectra::SparseGenMatProd<double>> eigs(op, 3, 1000);

    // Initialize and compute
    eigs.init();
    eigs.compute(Spectra::SortRule::LargestMagn);

    // Retrieve results
    Eigen::VectorXcd evalues;
    if(eigs.info() == Spectra::CompInfo::Successful) {
      evalues = eigs.eigenvalues();
      double lambda0 = std::abs(evalues(0));
      double lambda1 = std::abs(evalues(1));
      double lambda2 = std::abs(evalues(2));
      double gap     = lambda0 - lambda1;
      std::cout << log10(num_bits) << ","      
                << lambda0 << "," 
                << lambda1 << "," 
                << lambda2 << "," 
                << gap << "," 
                << log10(gap) << std::endl;
    } else {
      std::cout << "Spectral decomposition failed" << std::endl;
    }    

    num_bits *= 10;
    num_values /= 10;
  }

  return 0;    
}
