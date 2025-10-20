// This code is a companion to the paper "Negative Drift and State Instability 
// in a Bitwise System Equivalent to the Collatz Conjecture". This program 
// performs a statistical analysis of the bitwise dynamical system
// U(n) = 3n + nu(n), which is equivalent to the Collatz sequence.
// It runs a Monte Carlo simulation to:
// 1. Calculate the mean change in bit span (time average).
// 2. Generate transition matrices for a 16-state Markov model.
// 3. Validate the independence of MSB/LSB sub-systems via the Frobenius norm.
// 4. Calculate the stationary distribution and theoretical space averages.
//
// The Markov model states are defined by a 2-bit MSB state and 2-bit LSB state:
// LSB State: The 2 bits immediately to the left of the LSB
// MSB State: The 2 bits immediately to the right of the MSB
// State = 4 * (MSB State) + (LSB State)
//
// This code requires:
//   1. The GNU MP BigNum library (GMP), available here: https://gmplib.org/
//   2. The Eigen3 library for linear algebra, available here: https://eigen.tuxfamily.org
//
// To build:
// g++ -o collatz collatz.cpp -lgmp -O3
//
// You may need to specify the location of your Eigen3 headers:
// g++ -o collatz collatz.cpp -lgmp -O3 -I /path/to/eigen3
//
// Auth: Travis Galoppo (tgaloppo@gmail.com)
// Date: 2025-10-14
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <gmp.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/KroneckerProduct>

// minimum bits between LSB and MSB for
// inclusion in markov transition tables
#define MIN_TRANSITION_BITS 16

int main(int argc, char** argv) {
  mpz_t n, nu, two, tmp;
  gmp_randstate_t rnd_state;
  double Delta_L = 0.0;
  double Delta_H = 0.0;
  double Delta_S = 0.0;
  unsigned long long count = 0;
  bool use_transition;

  // default number of values and bit size
  size_t num_values = 1000000;
  size_t num_bits   = 128;

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
  size_t S0, S1;
  size_t msb0, lsb0;
  size_t state0, state1;

  // matrix buffers
  double Mlsb[4*4];
  double Mmsb[4*4];
  double Mfull[16*16];
  memset(Mlsb, 0, sizeof(double) * 4 * 4);
  memset(Mmsb, 0, sizeof(double) * 4 * 4);
  memset(Mfull, 0, sizeof(double) * 16 * 16);

  // main loop
  for (size_t j = 0; j < num_values; j++) {
    // show progress
    if (j % 1000 == 0) {
      if (j % 10000 == 0)
        std::cout << j;
      else
        std::cout << ".";
    } 
    std::cout << std::flush;   

    // generate a random starting value
    mpz_urandomb(n, rnd_state, num_bits);
  
    // find LSB, MSB and initial span
    L0 = mpz_scan1(n, 0);
    H0 = mpz_sizeinbase(n, 2);
    S0 = H0 - L0;

    // loop until we reach a power of 2
    while ( mpz_popcount(n) > 1 ) {
      // 2-adic value of n
      mpz_pow_ui(nu, two, L0);

      // state transitions are only counted if
      // there are at least MIN_TRANSITION_BITS
      // between the LSB and MSB
      use_transition = H0 - L0 >= MIN_TRANSITION_BITS;

      // collect starting state information
      if (use_transition) {
        mpz_tdiv_q_2exp(tmp, n, L0+1);
        lsb0 = mpz_get_ui(tmp) & 3; // 2 bits "left of" the LSB
        mpz_tdiv_q_2exp(tmp, n, H0-3);
        msb0 = mpz_get_ui(tmp) & 3; // 2 bits "right of" the MSB
        state0 = msb0 * 4 + lsb0;
      }

      // apply U(n)
      mpz_mul_ui(n, n, 3);
      mpz_add(n, n, nu);
    
      // find LSB, MSB, and span of resulting number
      L1 = mpz_scan1(n, L0); // L1 must be >L0, so this safe
      H1 = mpz_sizeinbase(n,2);
      S1 = H1 - L1;

      use_transition = use_transition && H1 - L1 >= MIN_TRANSITION_BITS;

      // (if applicable) collect transitioned state
      // info and record transition
      if (use_transition) {
        mpz_tdiv_q_2exp(tmp, n, L1+1);
        size_t lsb1 = mpz_get_ui(tmp) & 3;
        mpz_tdiv_q_2exp(tmp, n, H1-3);
        size_t msb1 = mpz_get_ui(tmp) & 3;
        Mlsb[lsb0 * 4 + lsb1] += 1.0;
        Mmsb[msb0 * 4 + msb1] += 1.0;
        state1 = msb1 * 4 + lsb1;
        Mfull[state0 * 16 + state1] += 1.0;
      }

      // Record changes in LSB, MSB, and span
      // (note this is not conditional on span size)
      Delta_L += L1 - L0;
      Delta_H += H1 - H0;
      Delta_S += (double)S1 - (double)S0;

      // Count total number of transitions
      count += 1;
   
      // Update starting conditions for next iteration
      L0 = L1;
      H0 = H1;
      S0 = S1;

      // if we want to histogram the fractional positions
      // generated, we can use this to get the values.
      /*
      if (mpz_popcount(n) > 1) {
        mpq_t a;
        mpq_init(a);
        mpq_set_z(a, n);
        mpq_div_2exp(a, a, H0);
        double alpha = mpq_get_d(a);
        std::cerr << alpha << std::endl;
        mpq_clear(a);
      }
      */
    }
  }
  std::cout << std::endl << std::endl;
  
  // Compute mean observed deltas
  double mean_delta_l = Delta_L / count;
  double mean_delta_h = Delta_H / count;
  double mean_delta_s = Delta_S / count;
  
  std::cout << "Total starting values: " << num_values << std::endl;
  std::cout << "         Max Bit size: " << num_bits << std::endl;
  std::cout << "    Total Transitions: " << count << std::endl;
  std::cout << std::setprecision(8);
  std::cout << "         Mean Delta_L: " << mean_delta_l << std::endl 
            << "         Mean Delta_H: " << mean_delta_h << std::endl
            << "         Mean Delta_S: " << mean_delta_s << std::endl 
            << std::endl;

  mpz_clears(n, nu, two, tmp, NULL);

  // Convert LSB transition counts to probabilities
  // (normalize rows to sum to 1)
  Eigen::Map<Eigen::Matrix<double, 4, 4>, Eigen::Unaligned, Eigen::Stride<1,4>> M_lsb(Mlsb, 4, 4);
  M_lsb = M_lsb.array().colwise() / M_lsb.array().rowwise().sum();

  // Convert MSB transition counts to probabilities
  Eigen::Map<Eigen::Matrix<double, 4, 4>, Eigen::Unaligned, Eigen::Stride<1,4>> M_msb(Mmsb, 4, 4);
  M_msb = M_msb.array().colwise() / M_msb.array().rowwise().sum();

  // Convert joint transition counts to probabilities
  Eigen::Map<Eigen::Matrix<double, 16, 16>, Eigen::Unaligned, Eigen::Stride<1,16>> M_full(Mfull, 16, 16);
  M_full = M_full.array().colwise() / M_full.array().rowwise().sum();

  // Display LSB and MSB transition tables
  std::cout << "LSB Transition Matrix:" << std::endl << M_lsb << std::endl << std::endl;
  std::cout << "MSB Transition Matrix:" << std::endl << M_msb << std::endl << std::endl;

  // Compute joint transition matrix from individual system matrices
  Eigen::Matrix<double, 16, 16> M_prime = Eigen::KroneckerProduct(M_msb, M_lsb);

  // Compute the Frobenius norm between the observed full transition
  // matrix and that computed from individual systems; if the systems
  // are independent, this should be very nearly zero.
  double norm = (M_prime - M_full).norm();
  std::cout << "Frobenius norm between M_full and M_prime: " << norm << std::endl << std::endl;

  // Identity matrix for stationary vector calculations
  Eigen::Matrix<double, 4, 4> Eye = Eigen::Matrix<double, 4, 4>::Identity();

  // Compute stationary vector of LSB system
  Eigen::FullPivLU<Eigen::Matrix<double, 4, 4>> lu_lsb(Eye - M_lsb.transpose());
  Eigen::Vector<double, 4> pi_lsb = lu_lsb.kernel();
  pi_lsb = pi_lsb / pi_lsb.sum();

  // Compute stationary vector of MSB system
  Eigen::FullPivLU<Eigen::Matrix<double, 4, 4>> lu_msb(Eye - M_msb.transpose());
  Eigen::Vector<double, 4> pi_msb = lu_msb.kernel();
  pi_msb = pi_msb / pi_msb.sum();

  std::cout << "Stationary Vectors of Markov Model:" << std::endl;
  std::cout << "pi_lsb = " << pi_lsb.transpose() << std::endl;
  std::cout << "pi_msb = " << pi_msb.transpose() << std::endl << std::endl;

  // These are the expected bit shifts for each state; 3/4 of them
  // are deterministic; only LSB[2] and MSB[1] are not deterministic.
  // See propositions 7-9 and section 7.4 in the paper.
  Eigen::Vector<double, 4> lsb_shifts { 2.0, 1.0, 4.0, 1.0 };
  Eigen::Vector<double, 4> msb_shifts { 1.0, 5.0 / 3.0, 2.0, 2.0 };

  // Compute space averages of system
  double E_Delta_L = pi_lsb.dot(lsb_shifts);
  double E_Delta_H = pi_msb.dot(msb_shifts);
  double E_Delta_S = E_Delta_H - E_Delta_L;

  std::cout << "Space Averages of Markov Model:" << std::endl
            << "E[Delta_L] = " << E_Delta_L << std::endl
            << "E[Delta_H] = " << E_Delta_H << std::endl
            << "E[Delta_S] = " << E_Delta_S << std::endl << std::endl;

  return 0;
}

