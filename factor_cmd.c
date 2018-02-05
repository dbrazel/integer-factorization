#include "rho.c"

int main() {
	mpz_t n, p, q, c;
	mpz_inits(n, p, q, c, NULL);
	int flag;
	
	printf("Enter a number to be factored: ");
	flag = mpz_inp_str(n, NULL, 10);
	assert(flag != 0);
	
	mpz_set_ui(c, 1);
	
	flag = rho(n, p, q, c);
	
	if (flag != 0) {
		printf("\nThe factorization failed.\n");
		return 1;
	}
	
	printf("The factors are \n");
	mpz_out_str(NULL, 10, p);
	printf("\nand \n");
	mpz_out_str(NULL, 10, q);
	printf("\n");
	
	mpz_clears(n, p, q, c, NULL);
	
	return 0;
}
