// This code is a companion to the paper "Deterministic Limits and Ergodic
// Properties of a Bitwise Syracuse Map". This program generates a frequency
// table for the LSB shifts observed in Syracuse trajectories.
//
// This code requires:
//   1. The GNU MP BigNum library (GMP), available here: https://gmplib.org/
//
// To build:
// g++ -o delta-L-freq delta-L-freq.cpp -lgmp -O3
//
// Auth: Travis Galoppo (tgaloppo@gmail.com)
// Date: 2025-11-19
#include <iostream>
#include <map>
#include <gmp.h>

int main(int argc, char** argv) {
  mpz_t n, nu, two, tmp;
  gmp_randstate_t rnd_state;

  // default number of values and bit size
  size_t num_values = 1;
  size_t num_bits   = 1000000;

  // command line args <num_values> <num_bits>
  if (argc > 1)
    num_values = std::atol(argv[1]);
  if (argc > 2)
    num_bits   = std::atol(argv[2]);

  // initialize
  gmp_randinit_default(rnd_state);
  gmp_randseed_ui(rnd_state, time(NULL));
  
  mpz_inits(n, nu, two, tmp, NULL);
  
  mpz_set_ui(two, 2);

  size_t L0, L1;
  size_t H0, H1;

  std::map<long, long long> freq_map;

  // main loop
  for (size_t j = 0; j < num_values; j++) {
    // generate a random starting value
    mpz_urandomb(n, rnd_state, num_bits);
  
    // find LSB, MSB and initial span
    L0 = mpz_scan1(n, 0);
    H0 = mpz_sizeinbase(n, 2) - 1;

    // loop until we reach a power of 2
    while ( H0 > L0 ) {
      // 2-adic value of n
      mpz_pow_ui(nu, two, L0);

      // apply U(n)
      mpz_mul_ui(n, n, 3);
      mpz_add(n, n, nu);
    
      // find LSB, MSB, and span of resulting number
      L1 = mpz_scan1(n, L0); // L1 must be >L0, so this safe
      H1 = mpz_sizeinbase(n,2) - 1;

      long dL = L1 - L0; 
      // Record it
      freq_map[dL]++;

      // Update starting conditions for next iteration;
      // for efficiency, pull LSB back to position 0
      mpz_tdiv_q_2exp(n, n, L1);
      L0 = 0;
      H0 = H1 - L1;
    }
  }

  mpz_clears(n, nu, two, tmp, NULL);

  std::cout << "dL, Count" << std::endl;
  for (auto const& [shift, count] : freq_map) {
    std::cout << shift << ", " << count << std::endl;
  }

  return 0;
}

