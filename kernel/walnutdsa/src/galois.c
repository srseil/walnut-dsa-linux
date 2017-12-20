#include <linux/log2.h>

#include "galois.h"

uint8_t galois_mult(unsigned int a, unsigned int b, unsigned int q);
uint8_t galois_inverse(unsigned int x, unsigned int q);

/* Irreducable polynomials for the finite fields GF(2^n) with n = 5, 6, 7, 8. */
static unsigned int irr_polys[4] = {
	37,  /* x^5 + x^2 + 1 */
	67,  /* x^6 + x + 1 */
	131, /* x^7 + x + 1 */
	283, /* x^8 + x^4 + x^3 + x + 1 */
};

uint8_t galois_mult(unsigned int a, unsigned int b, unsigned int q) {
	uint8_t result = 0;
	unsigned int poly = irr_polys[(unsigned int) ilog2(q) - 5];
	bool overflow;
	while (a != 0) {
		if ((a & 1) != 0)
			result ^= b;
		overflow = (b & (q / 2)) != 0;
		b <<= 1;
		if (overflow)
			b ^= poly;
		a >>= 1;
	}
	return result;
}


uint8_t galois_inverse(unsigned int x, unsigned int q) {
	int i;
	for (i = 0; i < q; i++) {
		if (galois_mult(x, i, q) == 1)
			return i;
	}
	return 0;
}

