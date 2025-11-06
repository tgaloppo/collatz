#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <gmp.h>
#include <eigen3/Eigen/Dense>
#include <cmath>

int main(int argc, char** argv) {
    const size_t n_orbits = 10000;
    size_t n_bits = 128;
    size_t L0;
    
    mpz_t n, nu, two, mask, r, n_shifted;
    gmp_randstate_t rnd_state;

    // look at the lowest k_bits bits
    const size_t k_bits = 16;
    // We'll divide by 2^k later. Pre-calculate this as a double.
    const double divisor_d = exp2(k_bits); // This is 2^16 = 65536.0

    mpz_inits(n, nu, two, mask, r, n_shifted, NULL);
    mpz_set_ui(two, 2);

    // Create the bitwise mask
    // mask = 2^k - 1
    mpz_set_ui(mask, 1);
    mpz_mul_2exp(mask, mask, k_bits); // mask = 2^k
    mpz_sub_ui(mask, mask, 1);      // mask = 2^k - 1 (e.g., 0b1111111111111111)

    gmp_randinit_default(rnd_state);
    gmp_randseed_ui(rnd_state, time(NULL));

    std::cout << "nbits,orbit,lag,autocorrelation" << std::endl;

    do {
        for (size_t orbit = 0; orbit < n_orbits; orbit++) {
            mpz_urandomb(n, rnd_state, n_bits);
            L0 = mpz_scan1(n, 0);

            std::vector<double> low_bits_frac;
            while ( mpz_popcount(n) > 1) {
                // n_shifted = n / 2^L0 (This is the integer right-shift to get the odd part)
                mpz_fdiv_q_2exp(n_shifted, n, L0);
                
                // r = n_shifted & (2^k - 1) (get the k lowest bits of the odd part)
                mpz_and(r, n_shifted, mask);
                
                // mpz_get_d is now safe, as 'r' is a small integer (0 to 2^k-1)
                // This gives the fractional value of the k lowest bits
                low_bits_frac.push_back( mpz_get_d(r) / divisor_d );

                // Collatz iteration
                mpz_pow_ui(nu, two, L0);
                mpz_mul_ui(n, n, 3);
                mpz_add(n, n, nu);
                L0 = mpz_scan1(n, L0);
            }

            if (low_bits_frac.size() < 50) continue; // Skip orbits that are too short for ACF

            double mu = 0.0;
            for (size_t j = 0; j < low_bits_frac.size(); j++) {
                mu += low_bits_frac[j];
            }
            mu /= low_bits_frac.size();

            double variance = 0.0;
            for (size_t j=0; j < low_bits_frac.size(); j++) {
                variance += (low_bits_frac[j] - mu) * (low_bits_frac[j] - mu);
            }
            variance /= low_bits_frac.size();

            if (variance == 0) continue; // Avoid division by zero

            // Cap max_lag to a reasonable number (e.g., 3000) or half the size
            size_t max_lag = low_bits_frac.size(); //std::min((size_t)3000, low_bits_frac.size() / 2);
            
            for (size_t lag = 1; lag < max_lag; lag++) {
                double sum = 0.0;
                for (size_t j=0; j < low_bits_frac.size() - lag; j++) {
                    sum += (low_bits_frac[j] - mu) * (low_bits_frac[j+lag] - mu);
                }
                double autocovariance  = sum / (low_bits_frac.size() - lag);
                double autocorrelation = autocovariance / variance;
                std::cout << n_bits << "," << orbit << "," << lag << "," << autocorrelation << std::endl;
            }
        }

        n_bits *= 2;
    } while (n_bits < 4*1024);

    mpz_clears(n, nu, two, mask, r, n_shifted, NULL);
}

