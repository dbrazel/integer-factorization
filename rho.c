#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include <assert.h>
#include <time.h>
#include <stdbool.h>

/* An implementation of Pollard's Rho algorithm, for integer 
 * factorization. Requires three parameters, all mpz_t objects: 
 * n - the integer to be factored; p and q - these will be set to 
 * the factors, if they are found; c - the parameter in the 
 * pseudorandom function. The function will return zero if non-trivial 
 * prime factors are found. It will return one if an incorrect factor
 * is found*/

int rho (const mpz_t n, mpz_t p, mpz_t q, mpz_t c)
{
	/* Initialize and set variables used for cycle detection */
	/*unsigned int i = 1;
	unsigned int k = 2;*/
	
	mpz_t i, k;
	mpz_inits(i, k, NULL);
	mpz_set_ui(i, 1);
	mpz_set_ui(k, 2);
	
	/* Initialize and set x and y to a random integer
	 * between 0 and n-1 */
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));
	
	mpz_t x, y;
	mpz_init(x);
	mpz_init(y);
	mpz_urandomm(x, state, n);
	mpz_set(y, x);

	/* Initialize and set other mpz_t objects needed */
	mpz_t diff, d, remain;
	mpz_init(diff);
	mpz_init(d);
	mpz_init(remain);

	while (true)
	{
		/*i++;*/
		mpz_add_ui(i, i, 1);
		mpz_pow_ui(x, x, 2);
		mpz_sub(x, x, c);
		mpz_mod(x, x, n);
		mpz_sub(diff, x, y);
		mpz_abs(diff, diff);
		mpz_gcd(d, diff, n);
		
		/*if (i == k)*/
		if (mpz_cmp(i, k) == 0)
		{
			mpz_set(y, x);
			/*k = 2*k;*/
			mpz_mul_ui(k, k, 2);
			/*printf("Iteration ");
			mpz_out_str(NULL, 10, i);
			printf("\n");*/
		}

		/* Check if d is equal to n, if so we have failed. Need to 
		 * rerun with a new c and new random start point*/
		if (mpz_cmp(d, n) == 0)
		{
			mpz_add_ui(c, c, 2);
			mpz_clears(x, y, diff, d, remain, i, k, NULL);
			printf("Rerunning...\n");
			return rho(n, p, q, c);
		}

		if (!(mpz_cmp_ui(d, 1) == 0))
		{
			/* We have a putative factor. First, check it*/
			mpz_tdiv_qr(q, remain, n, d);
			
			/* If we didn't get a factor, return with an error code */
			if (!(mpz_cmp_ui(remain, 0) == 0)) {
				mpz_clears(x, y, diff, d, remain, i, k, NULL);
				return 1;
			}
			
			mpz_set(p, d);
			mpz_clears(x, y, diff, d, remain, i, k, NULL);
			return 0;
		}
	}	
}
