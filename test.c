#define _GNU_SOURCE

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>

#include <sys/syscall.h>
#include <linux/random.h>

#include "galois.h"

#define PURE_BRAID_SIZE(low, high) (2 * ((high) - (low)) + 2)

/* Note that both the artin generators and the numbers in the permutations
 * start at 1 and not at 0!
 */

struct cb_pair {
	int n;				// Number of strands of the braid group dealt with.
	int q;				// Size of the galois field operated in.
	int *matrix;		// Pointer to an n*n matrix over the galois field.
	int *permutation;	// Pointer to an array of numbers modulo n.
};


void emult(struct cb_pair *result, struct cb_pair m, int braid, int *t_values);
void gen_cloaking_elem(int **dest, int *dest_size,
	int a, int b, int c, int n, int bits);
int *gen_priv_key(int a, int b, int n, int bits);
void gen_pure_braid(int *dest, int x, int y, int n);
void encode(int **dest, int *dest_size, char *message, int size, int n);

void example();
void cloaking_example();


int main(int argc, char *argv[]) {
	int *cloak = NULL;
	int cloak_size = 0;
	gen_cloaking_elem(&cloak, &cloak_size, 2, 5, 5, 8, 32);
	return 0;

	//cloaking_example();
	//int *priv = gen_priv_key(3, 6, 8, 128);
	//example();

	char msg[] = {121u, 160u}; // 0111 1001 1010 0000
	int *enc, enc_size;
	encode(&enc, &enc_size, msg, 2, 8);

	for (int i = 0; i < enc_size; i++)
		printf("%i, ", enc[i]);
	printf("\n");

	return 0;
}

void encode(int **dest, int *dest_size, char *message, int size, int n) {
	/* "Choose and fix a subset of four generators"
	 * Apparently I can pick arbitrary ones... maybe choose randomly?
	 */

	// Reserve enough space for the unreduced worst case.
	int *braid = calloc(PURE_BRAID_SIZE(1, n - 1) * 4 * size * 2, sizeof(char));
	unsigned char digest = *message;

	int mark = 0;            // Marker for the braid array index.
	int gen1, gen2,          // Generators for the digests, two low bits.
		power1, power2;      // Powers of the generators, two high bits.
	int reduced_middle = 0;  // Flag: Remove one reduced middle generator?
	int digest_complete = 0; // Flag: Get new digest from message?

	// Start with gen2 for better loop management.
	gen2 = ((48 & digest) >> 4) + 1;
	power2 = ((192 & digest) >> 6) + 1;
	for (int i = n - 1; i > gen2; i--, mark++) {
		braid[mark] = i;
	}
	printf("Initial: gen: %i, powert: %i\n", gen2, power2);

	// Iterate over the message in 4-bit increments.
	for (int i = 0; i < size * 2 - 1; i++) {
		// Update generators and powers.
		gen1 = gen2;
		power1 = power2;
		if (digest_complete) {
			digest = message[(i + 1) / 2];
			digest_complete = 0;
		} else {
			digest <<= 4;
			digest_complete = 1;
		}
		gen2 = ((48 & digest) >> 4) + 1;
		power2 = ((192 & digest) >> 6) + 1;

		printf("digest: %u\n", digest);
		printf("gen: %i, power: %i\n", gen2, power2);

		// Build first generator middle part.
		for (int i = 0 + reduced_middle; i < power1 * 2; i++, mark++) {
			braid[mark] = gen1;
		}

		// Build freely reduced start and end connecting the two generators.
		if (gen1 > gen2) {
			for (int i = gen1; i > gen2; i--, mark++) {
				braid[mark] = i;
			}
			reduced_middle = 0;
		} else if (gen1 < gen2) {
			for (int i = gen1 + 1; i < gen2; i++, mark++) {
				braid[mark] = i * -1;
			}
			reduced_middle = 1;
		}
	}
	
	// Build last generator middle part.
	for (int i = 0 + reduced_middle; i < power2 * 2; i++, mark++) {
		braid[mark] = gen2;
	}

	// Build last generator ending.
	for (int i = gen2 + 1; i <= n - 1; i++, mark++) {
		braid[mark] = i * -1;
	}

	*dest = braid;
	*dest_size = mark;
}

void gen_pure_braid(int *dest, int low, int high, int n) {
	/* Generates high, high-1, ..., low, low, ..., (high-1)^-1, high^-1 */
	int size = PURE_BRAID_SIZE(low, high);

	dest[high - low] = low;
	dest[high - low + 1] = low;
	for (int i = 0; high > low; i++, high--) {
		dest[i] = high;
		dest[size - 1 - i] = -1 * high;
	}

	for (int i = 0; i < size; i++)
		printf("%i, ", dest[i]);
	printf("\n");
}

int *gen_priv_key(int a, int b, int n, int bits) {
	int l = bits / (2 * log2f((n - 1) * (n - 2))) + 1;
	printf("Security level: %i, L = %i\n", bits, l);

	// Pick random delimiters x and y so that 1 <= x < y < n.
	int xs[2 * l], ys[2 * l];
	for (int i = 0; i < 2 * l; i++) {
		int rand = 0;

		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		rand = abs(rand);
		int x = rand % (n - 2) + 1;
		xs[i] = x;

		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		rand = abs(rand);
		int y = x + (rand % (n - x - 1)) + 1;
		ys[i] = y;
	}

	for (int i = 0; i < 2 * l; i++)
		printf("x: %i, y: %i\n", xs[i], ys[i]);
	printf("\n\n\n");

	// Calculate braid size
	int size = 1;
	for (int i = 0; i <  2 * l; i += 1)
		size += 2 * (ys[i] - xs[i]) + 2;

	printf("Private key braid size: %i\n", size);

	// Build braid
	int *braid = malloc(size * sizeof(int));
	int mark = 0;
	for (int i = 0; i < 2 * l; i++) {
		if (i == l) {
			int rand = 0;
			int middle = a;
			while (middle == a || middle == a - 1
					|| middle == b || middle == b - 1) {
				syscall(SYS_getrandom, &rand, sizeof(int), 0);
				middle = abs(rand) % (n - 1) + 1;
			}
			braid[mark] = middle;
			mark++;

			printf("middle: %i\n", middle);
		}

		gen_pure_braid(braid + mark, xs[i], ys[i], n);
		mark += 2 * (ys[i] - xs[i]) + 2;
	}

	printf("\n\n");

	for (int i = 0; i < size; i++)
		printf("%i, ", braid[i]);
	printf("\n");

	return braid;
}

void gen_cloaking_elem(int **dest, int *dest_size,
		int a, int b, int c, int n, int bits) {
	int mark = 0;
	int l = bits / (2 * log2f((n - 1) * (n - 2))) + 1;

	// Allocate enough space for the worst case.
	// This needs to be refined!
	int *braid = malloc(PURE_BRAID_SIZE(1, n - 1) * 4 * l * sizeof(int));


	// Build permuting part at start.
	if (c <= a) {
		for (int i = 0; i < a - c; i++, mark++)
			braid[mark] = a - i - 1;
		for (int i = 0; i < b - (c + 1); i++, mark++)
			braid[mark] = b - i - 1;
	} else {
		for (int i = 0; i < (c + 1) - b; i++, mark++)
			braid[mark] = b + i;
		for (int i = 0; i < c - a; i++, mark++)
			braid[mark] = a + i;
	}

	// Randomly build l pure braids.
	for (int i = 0; i < l; i++) {
		int rand = 0, high, low;
		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		high = abs(rand) % (n - 2) + 2; // 1 < high < 8
		do {
			syscall(SYS_getrandom, &rand, sizeof(int), 0);
			low = abs(rand) % (high) + 1; // 1 <= low < high
		} while (low == high);
		printf("Low: %i, High: %i\n", low, high);

		gen_pure_braid(&(braid[mark]), low, high, n);
		mark += PURE_BRAID_SIZE(low, high);
	}

	// Build middle part.
	braid[mark] = b;
	braid[mark + 1] = b;
	mark += 2;

	// Build inverse of the first half after the middle.
	int half = mark - 2;
	for (int i = 0; i < half; i++, mark++) {
		braid[mark] = -1 * braid[half - 1 - i];
	}

	*dest = braid;
	*dest_size = mark;

	for (int i = 0; i < mark; i++)
		printf("%i, ", braid[i]);
	printf("\n");
}

void emult(struct cb_pair *result, struct cb_pair m, int braid, int *t_values) {
	int n = m.n;		// Number of braid strands
	int q = m.q;		// Galois field size
	int w = log2(q);	// Prime potency for galois field arithmetic
	int abs_braid = abs(braid);

	printf("Artin generator: %i\n", braid);

	// T-variable permutation and evaluation
	int t_val;
	if (braid < 0) {
		int t_i = m.permutation[(abs_braid - 1 + 1) % n];
		t_val = galois_inverse(t_values[t_i - 1], w);
	} else {
		int t_i = m.permutation[braid - 1];
		t_val = t_values[t_i - 1];
	}

	printf("Evaluated t-value: %i\n", t_val);

	// Colored Burau matrix construction
	int braid_matrix[n * n];
	for (int row = 0; row < n; row++) {
		if (row == abs_braid - 1) {
			for (int col = 0; col < abs_braid - 2; col++)
				braid_matrix[row * n + col] = 0;
			if (braid < 0) {
				braid_matrix[row * n + abs_braid - 2] = 1;
				braid_matrix[row * n + abs_braid - 1] = t_val;
				if (abs_braid != n)
					braid_matrix[row * n + abs_braid] = t_val;
			} else {
				braid_matrix[row * n + braid - 2] = t_val;
				braid_matrix[row * n + braid - 1] = t_val;
				if (braid != n)
					braid_matrix[row * n + braid] = 1;
			}
			for (int col = abs_braid + 1; col < n; col++)
				braid_matrix[row * n + col] = 0;
		} else {
			for (int col = 0; col < n; col++)
				braid_matrix[row * n + col] = (row == col) ? 1 : 0;
		}
	}

	printf("Constructed Colored Burau matrix:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%2i, ", braid_matrix[i * n + j]);
		}
		printf("\n");
	}

	// Matrix multiplication
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				result->matrix[i * n + j] ^= galois_single_multiply(
					m.matrix[i * n + k], braid_matrix[k * n + j], w);
			}
		}
	}

	printf("Multiplied matrix:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%2i, ", result->matrix[i * n + j]);
		}
		printf("\n");
	}

	// Permutation composition
	for (int i = 0; i < n; i++) {
		if (i == abs_braid - 1)
			result->permutation[i] = m.permutation[i + 1];
		else if (i == abs_braid + 1 - 1)
			result->permutation[i] = m.permutation[i - 1];
		else
			result->permutation[i] = m.permutation[i];
	}

	printf("Composed permutation:\n");
	for (int i = 0; i < n; i++) {
		printf("%i -> %i,\n", i + 1, result->permutation[i]);
	}

	printf("\n");
}

void cloaking_example() {
	int *cloak, cloak_size;
	gen_cloaking_elem(&cloak, &cloak_size, 3, 6, 5, 8, 32);
	printf("Cloaking element size: %i\n", cloak_size);

	struct cb_pair cb;
	cb.n = 8;
	cb.q = 32;
	int cb_matrix[8*8] = {
		15, 30, 7, 18, 13, 20, 15, 31,
		10, 19, 13, 19, 6, 17, 11, 21,
		10, 14, 16, 19, 6, 17, 11, 21,
		17, 31, 21, 28, 15, 23, 2, 16,
		9, 16, 10, 13, 12, 7, 31, 20,
		9, 16, 10, 13, 21, 30, 31, 20,
		0, 0, 2, 4, 23, 0, 8, 24,
		0, 0, 0, 0, 0, 0, 0, 1,
	};
	cb.matrix = cb_matrix;
	int cb_permutation[8] = {2, 4, 3, 8, 1, 6, 7, 5};
	cb.permutation = cb_permutation;

	int t_values[8] = {28, 24, 1, 9, 26, 1, 18, 18};

	for (int i = 0; i < cloak_size; i++) {
		struct cb_pair result;
		result.n = 8;
		result.q = 32;
		int result_matrix[8*8] = {0};
		result.matrix = result_matrix;
		int result_permutation[8] = {0};
		result.permutation = result_permutation;

		emult(&result, cb, cloak[i], t_values);

		memcpy(cb.matrix, result.matrix, 8 * 8 * sizeof(int));
		memcpy(cb.permutation, result.permutation, 8 * sizeof(int));
	}
}

void example() {
	struct cb_pair m;
	m.n = 8;
	m.q = 32;

	int m_matrix[8*8] = {
		1, 0, 0, 0, 0, 0, 0, 0,
		0, 1, 0, 0, 0, 0, 0, 0,
		0, 0, 1, 0, 0, 0, 0, 0,
		0, 0, 0, 1, 0, 0, 0, 0,
		0, 0, 0, 0, 1, 0, 0, 0,
		0, 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0, 0, 0, 1
	};
	m.matrix = malloc(8 * 8 * sizeof(int));
	memcpy(m.matrix, &m_matrix, 8 * 8 * sizeof(int));

	int m_permutation[8] = {1, 2, 3, 4, 5, 6, 7, 8};
	m.permutation = malloc(8 * sizeof(int));
	memcpy(m.permutation, &m_permutation, 8 * sizeof(int));

	int braid1[] = {
		-4,3,-6,1,-4,-3,1,-4,-5,2,-1,-2,-2,-3,-1,5,5,-1,4,4,5,4,3,-4,-7,-3,6,5,
		-7,4,2,1,-7,-5,1,2,6,-4,-1,-2,-5,-1,4,-3,6,-3,5,1,-5,-2,4,-6,-7,-1,-1,3,
		-7,-4,-3,2,-5,2,5,-1,6,4,2,3,4,-3,2,3,-4,5,6,5,-4,2,-1,6,-7,-6
	};

	int braid2[] = {
		7,6,5,5,5,4,3,3,3,3,3,3,3,3,3,3,3,2,1,1,1,1,1,1,1,1,-2,-3,-4,5,5,5,-6,7,
		7,6,5,4,3,3,3,3,-4,5,-6,7,7,7,7,7,7,6,5,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
		3,3,3,3,3,3,3,3,3,-4,5,5,5,-6,7,7,7,7,7,7,7,7,7,7,6,5,4,3,3,-4,5,5,4,3,
		2,1,1,1,1,1,1,-2,3,3,2,1,1,1,1,1,1,-2,3,3,3,3,3,3,3,3,3,3,3,-4,5,5,5,5,
		5,5,5,-6,7,7,7,7,7,7,6,5,5,5,5,5,5,5,5,5,4,3,3,3,3,-4,5,5,5,-6,7,7,7,7,
		7,7,7,7,7,7,7,7,6,5,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,1,1,-2,3,-4,
		-5,-6,7,7,7,7,7,7,6,5,4,3,2,1,1,-2,3,-4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
		-6,7,7,7,7,7,7,7,7,6,5,4,3,2,1,1,-2,-3,-4,-5,-6,7,7,6,5,5,5,5,5,5,5,5,5,
		4,3,3,3,3,3,3,3,3,3,2,1,1,1,1,1,1,1,1,-2,-3,-4,5,5,5,5,5,5,5,5,5,5,4,3,
		3,3,3,3,3,3,3,-4,5,5,5,5,5,5,5,-6,7,7,7,7,7,7,7,7,6,5,4,3,2,1,1,1,1,1,1,
		1,1,1,1,1,1,-2,-3,-4,5,5,5,5,5,5,5,5,4,3,2,1,1,-2,-3,-4,5,5,5,5,4,3,2,1,
		1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-2,-3,-4,-5,-6,-7
	};

	int braid[] = {
		5, 4, 3, 3, -4, -5, 6, 5, 5, -6, 7, 6, 6, -7, 5, 4, 4, -5, 3, 2, 2, -3,
		6, 5, 5, -6, 7, 6, 5, 5, -6, -7, 7, 6, 6, -7, 7, 6, 6, -7, 4, 3, 3, -4,
		6, 5, 5, -6, 7, 6, 5, 4, 3, 2, 2, -3, -4, -5, -6, -7, 4, 3, 2, 2, -3, 6,
		5, 4, 3, 3, -4, -5, -6, 6, 5, 4, 4, -5, -6, 7, 6, 6, -7, 4, 3, 3, -4, 4,
		3, 2, 1, 1, -2, -3, -4, 3, 2, 1, 1, -2, -3, 6, 5, 5, -6, 6, 5, 5, -6, 5,
		4, 3, 2, 1, 1, -2, -3, -4, -5, 7, 6, 6, -7, 7, 6, 6, -7	
	};

	int tvals[8] = {28, 24, 1, 9, 26, 1, 18, 18};

	for (int i = 0; i < sizeof(braid) / sizeof(int); i++) {
		struct cb_pair result;
		result.n = 8;
		result.q = 32;
		int r_matrix[8*8] = {0};
		result.matrix = r_matrix;
		int r_permutation[8] = {0};
		result.permutation = r_permutation;

		emult(&result, m, braid[i], tvals);

		memcpy(m.matrix, result.matrix, 8 * 8 * sizeof(int));
		memcpy(m.permutation, result.permutation, 8 * sizeof(int));
	}

	free(m.matrix);
	free(m.permutation);
}

/*
void example2() {
	struct cb_pair2 m;
	m.n = 4;
	m.q = 32;
	int m_matrix[4*4] = {
		1, 2, 4, 1, 
		0, 1, 2, 0, 
		0, 3, 1, 0, 
		0, 0, 0, 1
	};
	m.matrix = malloc(4 * 4 * sizeof(int));
	memcpy(m.matrix, &m_matrix, 4 * 4 * sizeof(int));
	int m_permutation[4] = {2, 3, 1, 4};
	m.permutation = malloc(4 * sizeof(int));
	memcpy(m.permutation, &m_permutation, 4 * sizeof(int));

	int braid[] = {
		3, 3, 2, 1, -3, -2, 1, 1
	};

	int tvals[4] = {3, 2, 3, 4};
	for (int i = 0; i < sizeof(braid) / sizeof(int); i++) {
		struct m_pair result;
		result.n = 4;
		result.q = 5;
		int r_matrix[4*4] = {0};
		result.matrix = r_matrix;
		int r_permutation[4] = {0};
		result.permutation = r_permutation;

		emult(&result, m, braid[i], tvals);

		memcpy(m.matrix, result.matrix, 4 * 4 * sizeof(int));
		memcpy(m.permutation, result.permutation, 4 * sizeof(int));
	}

	free(m.matrix);
	free(m.permutation);
}
*/
