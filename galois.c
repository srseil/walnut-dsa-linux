#ifndef GALOIS
#define GALOIS

uint8_t galois_mult(uint8_t a, uint8_t b);
uint8_t galois_inverse(uint8_t x);

uint8_t galois_mult(uint8_t a, uint8_t b) {
	uint8_t r = 0, t;
	while (a != 0) {
		if ((a & 1) != 0)
			r ^= b;
		t = b & 0x10;
		b <<= 1;
		if (t != 0)
			b ^= 0x25;
		a >>= 1;
	}
	return r;
}


uint8_t galois_inverse(uint8_t x) {
	for (int i = 0; i < 32; i++) {
		if (galois_mult(x, i) == 1u)
			return i;
	}
	return 0;
}

#endif

