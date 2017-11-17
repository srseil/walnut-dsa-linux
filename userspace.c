#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/syscall.h>

#include "galois.h"
#include "walnut.h"
#include "debug.h"

int verify_pub_params(struct pub_params *params);
void gen_priv_key(int8_t **dest, int *dest_size,
		struct pub_params *params, int sec_level);
void gen_pub_key(uint8_t **dest_matrix, uint8_t **dest_perm,
		int8_t *priv_key, int priv_key_len, struct pub_params *params);
static unsigned long long fac(int x);
static void emult(uint8_t *matrix, uint8_t *perm, int8_t artin,
		uint8_t *t_values, int q, int n);

int main()
{
	int8_t *priv_key;
	int priv_key_len;
	uint8_t *pub_key_matrix, *pub_key_perm;
	uint8_t *t_values;
	int8_t generators[4];
	struct pub_params pub_params;
	int sec_level;

	int one1, one2;
	int rand, i;

	int8_t *signature;
	int signature_len;
	uint8_t hash[] = {
		0x79, 0xb7, 0xac, 0x30, 0x39, 0x4a, 0xff, 0x82,
		0x92, 0xed, 0x52, 0xf3, 0x8d, 0x52, 0x0c, 0xf8,
		0x2b, 0x01, 0xc8, 0x00, 0x63, 0x07, 0x46, 0xef,
		0x61, 0x4b, 0x26, 0xf0, 0x45, 0xc1, 0x09, 0x5d
	};
	int hash_size = sizeof(hash) / sizeof(uint8_t);

	pub_params.braid_group = 8;
	pub_params.galois_order = 256;

	t_values = malloc(pub_params.braid_group * sizeof(uint8_t));
	for (i = 0; i < pub_params.braid_group; i++) {
		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		t_values[i] = abs(rand) % (pub_params.galois_order - 1) + 1;
	}
	syscall(SYS_getrandom, &rand, sizeof(int), 0);
	one1 = abs(rand) % (pub_params.braid_group - 2) + 1;
	t_values[one1] = 1;
	one2 = one1;
	do {
		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		one2 = abs(rand) % (pub_params.braid_group - 2) + 1;
	} while (one1 == one2);
	t_values[one2] = 1;
	pub_params.t_values = t_values;

	generators[0] = 5;
	generators[1] = 4;
	generators[2] = 7;
	generators[3] = 2;
	pub_params.generators = generators;

	/*
	if (verify_pub_params(&params) != 0) {
		return -1;
	}
	*/

	sec_level = 128;

	/*
	int8_t priv_inline[] = {4, 5, -2, 7, -5, 7, -4, 7, -6, -6, 7, 1, 1, 2, 4, -5, 7, 7, 6, 6, 1, -5, 7, 5, 3, 1, 5, 2, 5, 7, 5, -4, 5, -2, -4, 2, 4, -1, -1, -6, -4, -1, -2, 7, -5, 1, -2, 4, 6, -1, -5, 6, 2, -1, 6, -7, -2, -3, -2, -5, -4, -7, -1, 2, 2, 4, -6, 2, -4, 1, 1, 7, 6, -5, 3, 7, 2, -6, -5, 1, 5, 4, 3, 6, -7, -2, 4, 4, -6, 3, -5, 6, 3, 7, -5, -4, -7, -1, -6, -6, 7, -6, -5, -7, 5, 3, 5, 5, 1, 4, -6, -5, 4, 1};
	priv_key = priv_inline;
	priv_key_len = sizeof(priv_inline) / sizeof(int8_t);

	uint8_t t_inline[] = {28,24,1,9,26,1,18,18};
	pub_params.t_values = t_inline;
	*/

	gen_priv_key(&priv_key, &priv_key_len, &pub_params, sec_level);

	gen_pub_key(&pub_key_matrix, &pub_key_perm,
		priv_key, priv_key_len, &pub_params);

	// SIGNATURE GENERATION AND VERIFICATION
	generate_signature(&signature, &signature_len, hash, hash_size, priv_key,
		priv_key_len, pub_key_perm, &pub_params, sec_level);

	verify_signature(signature, signature_len, hash, hash_size,
			pub_key_matrix, pub_key_perm, &pub_params);

	free(signature);
	free(pub_key_perm);
	free(pub_key_matrix);
	free(priv_key);
	free(t_values);
	return 0;
}

/* Verifies a set of public parameters for their correctness. */
int verify_pub_params(struct pub_params *params)
{
	int ones = 0;
	int i, j;

	if (params->braid_group < 8 || params->braid_group > 128) {
		fprintf(stderr, "Braid group not supported. "
			"Valid values are between 8 and 128.\n");
		return -1;
	} else if (params->galois_order != 32 && params->galois_order != 64
			&& params->galois_order != 128 && params->galois_order != 256) {
		fprintf(stderr, "Galois field order not supported. "
			"Valid values are 32, 64, 128 and 256.\n");
		return -1;
	} else if (params->t_values == NULL) {
		fprintf(stderr, "No t-values given.\n");
		return -1;
	}

	for (i = 0; i < params->braid_group; i++) {
		if (params->t_values[i] == 0
				|| params->t_values[i] >= params->galois_order) {
			fprintf(stderr, "Invalid t-value '%i'. "
				"Valid values are between 1 and %i\n",
				i, params->galois_order - 1);
			return -1;
		}
		ones++;
	}
	if (ones < 2) {
		fprintf(stderr, "Not enough ones in t-values. "
			"There must be at least two ones among the values.\n");
		return -1;
	}

	for (i = 0; i < 4; i++) {
		if (params->generators[i] < 1
				|| params->generators[i] >= params->braid_group) {
			fprintf(stderr, "Invalid Artin generator '%i'. "
				"Valid values are between 1 and %i.\n",
				i, params->braid_group - 1);
			return -1;
		}
		for (j = 0; j < i; j++) {
			if (params->generators[j] == params->generators[i]) {
				fprintf(stderr, "Duplicate Artin generator '%i'.\n", i);
				return -1;
			}
		}
	}
	
	return 0;
}

/* Generates a private key of the suitable minimum length for the specified
 * security level in bits, compatible with the specified public parameters. */
void gen_priv_key(int8_t **dest, int *dest_size,
		struct pub_params *params, int sec_level)
{
	int n = params->braid_group;
	float x_prev, x_curr;
	int len;
	int8_t *braid;
	int rand;
	int i;

	// Get minimum key length with Newton's method.
	x_prev = sec_level;
	for (i = 0; i < 3; i++) {
		x_curr = x_prev
			- (x_prev + (n - 2) * log2(x_prev) - (sec_level + log2(fac(n - 1))))
			/ ((n - 2) / (x_prev * log(2)) + 1);
		x_prev = x_curr;
	}
	len = (int) x_curr + 1;

	braid = malloc(len * sizeof(int8_t));
	for (i = 0; i < len; i++) {
		do {
			syscall(SYS_getrandom, &rand, sizeof(int), 0);
			braid[i] = abs(rand) % (2 * n - 1) - (n - 1);
		} while (braid[i] == 0 || (i > 0 && braid[i - 1] * -1 == braid[i]));
	}

	*dest = braid;
	*dest_size = len;
}

/* Generates the according public key for a specified private key, consisting
 * of a matrix and a permutation part and compatible with the specified public
 * parameters. */
void gen_pub_key(uint8_t **dest_matrix, uint8_t **dest_perm,
		int8_t *priv_key, int priv_key_len, struct pub_params *params)
{
	int n = params->braid_group;
	uint8_t *matrix, *perm;
	int i;

	matrix = malloc(n * n * sizeof(uint8_t));
	for (i = 0; i < n * n; i++)
		matrix[i] = (i % (n + 1) == 0) ? 1 : 0;

	perm = malloc(n * sizeof(uint8_t));
	for (i = 0; i < n; i++)
		perm[i] = i + 1;

	for (i = 0; i < priv_key_len; i++) {
		emult(matrix, perm, priv_key[i],
			params->t_values, params->galois_order, params->braid_group);
	}

	*dest_matrix = matrix;
	*dest_perm = perm;
}

/* Calculates the factorial of an integer. */
unsigned long long fac(int x) {
	unsigned long long res = x;
	int i;
	for (i = x - 1; i > 0; i--)
		res *= i;
	return res;
}

static void emult(uint8_t *matrix, uint8_t *perm, int8_t artin,
		uint8_t *t_values, int q, int n)
{
	int abs_artin = abs(artin);
	uint8_t *result_matrix = calloc(n * n, sizeof(uint8_t));
	uint8_t *result_perm = calloc(n, sizeof(uint8_t));
	uint8_t artin_matrix[n * n];
	uint8_t t_val;
	int row, col;
	int i, j, k;

	printf("Artin generator: %i\n", artin);

	// T-variable permutation and evaluation
	if (artin < 0) {
		int t_i = perm[(abs_artin - 1 + 1) % n];
		t_val = galois_inverse(t_values[t_i - 1], q);
	} else {
		int t_i = perm[artin - 1];
		t_val = t_values[t_i - 1];
	}
	printf("Evaluated t-value: %i\n", t_val);

	// Colored Burau matrix construction
	for (row = 0; row < n; row++) {
		if (row == abs_artin - 1) {
			for (col = 0; col < abs_artin - 2; col++)
				artin_matrix[row * n + col] = 0;
			if (artin < 0) {
				artin_matrix[row * n + abs_artin - 2] = 1;
				artin_matrix[row * n + abs_artin - 1] = t_val;
				if (abs_artin != n)
					artin_matrix[row * n + abs_artin] = t_val;
			} else {
				artin_matrix[row * n + artin - 2] = t_val;
				artin_matrix[row * n + artin - 1] = t_val;
				if (artin != n)
					artin_matrix[row * n + artin] = 1;
			}
			for (col = abs_artin + 1; col < n; col++)
				artin_matrix[row * n + col] = 0;
		} else {
			for (col = 0; col < n; col++)
				artin_matrix[row * n + col] = (row == col) ? 1 : 0;
		}
	}

	printf("Constructed Colored Burau matrix:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%2i, ", artin_matrix[i * n + j]);
		}
		printf("\n");
	}

	// Matrix multiplication
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				result_matrix[i * n + j] ^= galois_mult(
					matrix[i * n + k], artin_matrix[k * n + j], q);
			}
		}
	}

	printf("Multiplied matrix:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%2i, ", result_matrix[i * n + j]);
		printf("\n");
	}

	// Permutation composition
	for (i = 0; i < n; i++) {
		if (i == abs_artin - 1)
			result_perm[i] = perm[i + 1];
		else if (i == abs_artin + 1 - 1)
			result_perm[i] = perm[i - 1];
		else
			result_perm[i] = perm[i];
	}

	printf("Composed permutation:\n");
	for (i = 0; i < n; i++)
		printf("%i ", i + 1);
	printf("\n");
	for (i = 0; i < n; i++)
		printf("%i ", result_perm[i]);
	printf("\n\n");

	memcpy(matrix, result_matrix, n * n * sizeof(uint8_t));
	memcpy(perm, result_perm, n * sizeof(uint8_t));
	free(result_matrix);
	free(result_perm);
}
