#include <iostream>
#include <gmp.h>

int main(int argc, char** argv) {
    const size_t n_orbits = 2000;
    size_t n_bits = 32;
    size_t L0; 
    mpz_t n, nu, two;
    mpq_t a;
    gmp_randstate_t rnd_state;
    
    mpz_inits(n, nu, two, NULL);
    mpq_init(a);

    mpz_set_ui(two, 2);

    gmp_randinit_default(rnd_state);
    gmp_randseed_ui(rnd_state, time(NULL));

    std::cout << "nbits,orbit,alpha" << std::endl;
    do {
        for (size_t j = 0; j < n_orbits * (1024 / n_bits); j++) {
            // generate a random starting value
            mpz_urandomb(n, rnd_state, n_bits);  
            
            L0 = mpz_scan1(n, 0);

            while ( mpz_popcount(n) > 1 ) {
                // 2-adic value of n
                mpz_pow_ui(nu, two, L0); 

                // apply U(n)
                mpz_mul_ui(n, n, 3);
                mpz_add(n, n, nu);

                if (mpz_popcount(n) > 1) {
                    mpq_set_z(a, n);
                    mpq_div_2exp(a, a, mpz_sizeinbase(n,2));
                    double alpha = mpq_get_d(a);
                    std::cout << n_bits << "," << j << "," << alpha << std::endl;
                }  
                
                L0 = mpz_scan1(n, L0);
            }
        }

        n_bits *= 2;
    } while (n_bits < 2048);

    mpq_clear(a);
    mpz_clears(n, nu, two, NULL);

    return 0;
}
