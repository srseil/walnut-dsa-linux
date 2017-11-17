#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/syscall.h>

/*
#include <linux/kernel.h> for abs()
#include <linux/log2.h> for log2()
 */

#include "galois.h"
#include "bkl.h"
#include "dehornoy.h"
#include "debug.h"

#include "walnut.h"

#define PURE_BRAID_SIZE(x, n) (2 * ((n) - 1 - (x)) + 2)

void generate_signature(int8_t **dest, int *dest_len,
		uint8_t *hash, int hash_size,
		int8_t *priv_key, int priv_key_len,
		uint8_t *pub_key_perm,
		struct pub_params *params, int sec_level);
int verify_signature(int8_t *sig, int sig_len,
		uint8_t *hash, int hash_size,
		uint8_t *pub_key_matrix, uint8_t *pub_key_perm,
		struct pub_params *params);

//static int verify_pub_params(struct pub_params *params);
static void gen_pub_key(uint8_t **dest_matrix, uint8_t **dest_perm,
		int8_t *priv_key, int priv_key_len, struct pub_params *params);
static void encode(int8_t **dest, int *dest_size,
		uint8_t *msg, int msg_size,
		int8_t *gens, int n);
static void gen_pure_braids(int8_t **dest, int *dest_size, int min_len, int n);
static void gen_cloaking_elem(int8_t **dest, int *dest_size,
		uint8_t *perm, int a, int b, int n, int sec_level);
static void emult(uint8_t *matrix, uint8_t *perm, int8_t artin,
		uint8_t *t_values, int q, int n);

#if 0
int main()
{
#if 0
	//CREATE VALUE TABLE FOR GALOIS INVERSE IN G_256
	for (int i = 0; i < 256; i++) 
		printf("%i, ", galois_inverse(i));
	puts("\n");
	return 0;
#endif

#if 0
	// RANDOMIZED TEST FOR GALOIS ARITHMETIC IMPLEMENTATION
	while (true) {
		int rand;
		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		int x = abs(rand) % 32;
		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		int y = abs(rand) % 32;

		/*
		int a = galois_single_multiply(x, y, log2(32));
		//unsigned int b = galois_mult((unsigned int) x, (unsigned int) y);
		uint8_t b = galois((uint8_t) x, (uint8_t) y);
		printf("(%i * %i) -> %i == %u\n", x, y, a, b);
		*/

		int a = galois_inverse(x, log2(32));
		int b = glinverse(x);
		printf("%i == %i\n", a, b);
	}
	return 0;
#endif

#if 0
	// PICK FOUR RANDOM GENERATORS
	int duplicate;
	for (int i = 0; i < 4; i++) {
		do {
			int rand = 0;
			syscall(SYS_getrandom, &rand, sizeof(int), 0);
			gens[i] = abs(rand) % (n - 1) + 1;

			duplicate = 0;
			for (int j = 0; j < i; j++) {
				if (gens[j] == gens[i]) {
					duplicate = 1;
					break;
				}
			}
		} while (duplicate);
		printf("Randomly chosen generator: %i\n", gens[i]);
	}
#endif

#if 0
	int i;

	// Create public parameters struct.
	struct pub_params params;
	params.braid_group = 8;
	params.galois_order = 32;

	int rand = 0;
	uint8_t *t_values = malloc(params.braid_group * sizeof(uint8_t));
	for (i = 0; i < params.braid_group; i++) {
		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		t_values[i] = abs(rand) % (params.galois_order - 1) + 1;
	}
	syscall(SYS_getrandom, &rand, sizeof(int), 0);
	int one1 = abs(rand) % (params.braid_group - 2) + 1;
	t_values[one1] = 1;
	int one2 = one1;
	do {
		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		one2 = abs(rand) % (params.braid_group - 2) + 1;
	} while (one1 == one2);
	t_values[one2] = 1;

	params.t_values = t_values;
	int8_t generators[4] = {1, 3, 5, 7};
	params.generators = generators;
	if (verify_pub_params(&params) != 0) {
		return -1;
	}

	int sec_level = 128;

	int8_t *priv_key;
	int priv_key_len;
	gen_priv_key(&priv_key, &priv_key_len, &params, sec_level);

	uint8_t *pubk_matrix, *pubk_perm;
	gen_pub_key(
		&pubk_matrix, &pubk_perm, priv_key, priv_key_len, &params);

	int8_t *signature;
	int signature_len;

	uint8_t hash[] = {
		0x79, 0xb7, 0xac, 0x30, 0x39, 0x4a, 0xff, 0x82,
		0x92, 0xed, 0x52, 0xf3, 0x8d, 0x52, 0x0c, 0xf8,
		0x2b, 0x01, 0xc8, 0x00, 0x63, 0x07, 0x46, 0xef,
		0x61, 0x4b, 0x26, 0xf0, 0x45, 0xc1, 0x09, 0x5d
	};
	int hash_size = sizeof(hash) / sizeof(uint8_t);

	generate_signature(&signature, &signature_len, hash, hash_size, priv_key,
		priv_key_len, pubk_perm, &params, sec_level);

	verify_signature(signature, signature_len, hash, hash_size,
			pubk_matrix, pubk_perm, &params);

	free(signature);
	free(pubk_perm);
	free(pubk_matrix);
	free(priv_key);
	free(t_values);
#endif

	return 0;
}
#endif

#if 0
/* Verifies a set of public parameters for their correctness. */
int verify_pub_params(struct pub_params *params)
{
	int ones = 0;
	int i, j;

	if (params->braid_group < 8 || params->braid_group > 128) {
		// braid group not supported
		return -1;
	} else if (params->galois_order != 32 && params->galois_order != 64
			&& params->galois_order != 128 && params->galois_order != 256) {
		// galois field order not supported
		return -1;
	} else if (params->t_values == NULL) {
		// no t_values given
		return -1;
	}

	for (i = 0; i < params->braid_group; i++) {
		if (params->t_values[i] == 0
				|| params->t_values[i] >= params->galois_order) {
			// invalid t_value (not 0 because it needs to be invertable!)
			return -1;
		}
		ones++;
	}
	if (ones < 2) {
		// t_values are missing ones for a and b
		return -1;
	}

	for (i = 0; i < 4; i++) {
		if (params->generators[i] < 1
				|| params->generators[i] >= params->braid_group) {
			// invalid generator
			return -1;
		}
		for (j = 0; j < i; j++) {
			if (params->generators[j] == params->generators[i]) {
				// duplicate generator
				return -1;
			}
		}
	}
	
	return 0;
}
#endif

/* Generates the according public key for a specified private key, consisting
 * of a matrix and a permutation part and compatible with the specified public
 * parameters. */
static void gen_pub_key(uint8_t **dest_matrix, uint8_t **dest_perm,
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

void generate_signature(int8_t **dest, int *dest_len,
		uint8_t *hash, int hash_size,
		int8_t *priv_key, int priv_key_len,
		uint8_t *pub_key_perm,
		struct pub_params *params, int sec_level)
{
	/* Choose a and b randomly based on the occurencies of 1s in t_values? */
	int n = params->braid_group;

	int a = 0, b = 0;

	//int rand = 0; // Helper variable for random number generation.
	int mark = 0; // Helper variable for braid construction.

	int8_t *v, *v1, *v2; // Cloaking elements.
	int v_len, v1_len, v2_len; // Cloaking elements' lengths.

	uint8_t *identity_perm; // Identity permutation for pub key generation.

	int8_t *enc, *sig, *sig_old; // Encoded message and signature.
	int enc_len, sig_len; // Lengths of encoded message and signature.

	int i;

	identity_perm = malloc(n * sizeof(uint8_t));
	for (i = 0; i < n; i++)
		identity_perm[i] = i + 1;

	// Choose a and b based on t_values.
	for (i = 1; i < n - 1; i++) {
		if (params->t_values[i] == 1) {
			if (a == 0) {
				a = i + 1;
			} else {
				b = i + 1;
				break;
			}
		}
	}
	printf("A AND B: %i, %i", a, b);

	// Generate cloaking elements.
	gen_cloaking_elem(&v, &v_len, identity_perm, a, b, n, sec_level);
	gen_cloaking_elem(&v1, &v1_len, pub_key_perm, a, b, n, sec_level);
	gen_cloaking_elem(&v2, &v2_len, pub_key_perm, a, b, n, sec_level);
	free(identity_perm);

	// Encode hash into braid.
	encode(&enc, &enc_len, hash, hash_size, params->generators, n);
	
	// Build signature.
	sig_len = v_len +  2 * v1_len + v2_len + 2 * priv_key_len + enc_len;
	sig = malloc(sig_len * sizeof(int8_t));
	
	memcpy(sig, v2, v2_len * sizeof(int8_t));
	mark += v2_len;

	for (i = 0; i < v1_len; i++)
		sig[mark + i] = v1[v1_len - i - 1] * -1;
	mark += v1_len;

	for (i = 0; i < priv_key_len; i++)
		sig[mark + i] = priv_key[priv_key_len - i - 1] * -1;
	mark += priv_key_len;

	memcpy(sig + mark, v, v_len * sizeof(int8_t));
	mark += v_len;

	memcpy(sig + mark, enc, enc_len * sizeof(int8_t));
	mark += enc_len;

	memcpy(sig + mark, priv_key, priv_key_len * sizeof(int8_t));
	mark += priv_key_len;

	memcpy(sig + mark, v1, v1_len * sizeof(int8_t));
	mark += v1_len;

	free(enc);
	free(v2);
	free(v1);
	free(v);
	assert(mark == sig_len);

	// Convert signature into left canonical form.
	sig_old = sig;
	left_canonical_form(&sig, &sig_len, sig, sig_len, n);
	free(sig_old);

	// Apply Dehornoy's algorithm to signature.
	dehornoy(&sig, &sig_len, sig, sig_len);

	*dest = sig;
	*dest_len = sig_len;
}

int verify_signature(int8_t *sig, int sig_len,
		uint8_t *hash, int hash_size,
		uint8_t *pub_key_matrix, uint8_t *pub_key_perm,
		struct pub_params *params)
{
	int n = params->braid_group;
	int8_t *enc;
	int enc_len;
	uint8_t *msg_pub_key_matrix, *msg_pub_key_perm;
	uint8_t *matrix;
	int i, j, k;

	// Encode hash into braid.
	encode(&enc, &enc_len, hash, hash_size, params->generators, n);

	// Generate public key for encoded message.
	gen_pub_key(&msg_pub_key_matrix, &msg_pub_key_perm,
		enc, enc_len, params);
	free(enc);

	// Matrix multiplication
	matrix = calloc(n * n, sizeof(uint8_t));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				matrix[i * n + j] ^= galois_mult(msg_pub_key_matrix[i * n + k],
					pub_key_matrix[k * n + j], params->galois_order);
			}
		}
	}
	free(msg_pub_key_matrix);
	free(msg_pub_key_perm);

	for (i = 0; i < sig_len; i++) {
		emult(pub_key_matrix, pub_key_perm, sig[i], params->t_values,
			params->galois_order, params->braid_group);
	}

	printf("Multiplied matrix:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%2i, ", pub_key_matrix[i * n + j]);
		}
		printf("\n");
	}

	printf("Should be equal to:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%2i, ", matrix[i * n + j]);
		}
		printf("\n");
	}
	
	printf("\n");
	for (i = 0; i < n * n; i++) {
		if (matrix[i] != pub_key_matrix[i]) {
			free(matrix);
			printf("Verification FAILED! (at %i,%i [%i])\n", i%n, i/n, i);
			return -1;
		}
	}

	free(matrix);
	printf("Verification SUCCEEDED!\n");
	return 0;
}

/* Encodes a message given as a bit stream into a braid. */
static void encode(int8_t **dest, int *dest_size,
		uint8_t *msg, int msg_size,
		int8_t *gens, int n)
{
	// Reserve enough space for the unreduced worst case.
	int8_t *braid = calloc(PURE_BRAID_SIZE(1, n) * 4 * msg_size * 2,
		sizeof(int8_t));
	uint8_t digest = msg[0];
	int i, j;
	int mark = 0;            // Marker for the braid array index.
	int gen1, gen2,          // Generators for the digests, two low bits.
		power1, power2;      // Powers of the generators, two high bits.
	int reduced_middle = 0;  // Flag: Remove one reduced middle generator?
	int digest_complete = 0; // Flag: Get new digest from message?

	// Start with gen2 for better loop management.
	gen2 = gens[(48 & digest) >> 4];
	power2 = ((192 & digest) >> 6) + 1;
	for (i = n - 1; i > gen2; i--, mark++) {
		braid[mark] = i;
	}

	// Iterate over the message in 4-bit increments.
	for (i = 0; i < msg_size * 2 - 1; i++) {
		// Update generators and powers.
		gen1 = gen2;
		power1 = power2;
		if (digest_complete) {
			digest = msg[(i + 1) / 2];
			digest_complete = 0;
		} else {
			digest <<= 4;
			digest_complete = 1;
		}
		gen2 = gens[(48 & digest) >> 4];
		power2 = ((192 & digest) >> 6) + 1;

		// Build first generator middle part.
		for (j = 0 + reduced_middle; j < power1 * 2; j++, mark++)
			braid[mark] = gen1;

		// Build freely reduced start and end connecting the two generators.
		if (gen1 >= gen2) {
			for (j = gen1; j > gen2; j--, mark++)
				braid[mark] = j;
			reduced_middle = 0;
		} else if (gen1 < gen2) {
			for (j = gen1 + 1; j < gen2; j++, mark++)
				braid[mark] = j * -1;
			reduced_middle = 1;
		}
	}
	
	// Build last generator middle part.
	for (i = 0 + reduced_middle; i < power2 * 2; i++, mark++)
		braid[mark] = gen2;

	// Build last generator ending.
	for (i = gen2 + 1; i <= n - 1; i++, mark++)
		braid[mark] = i * -1;

	*dest = braid;
	*dest_size = mark;
}

static void gen_pure_braids(int8_t **dest, int *dest_size, int min_len, int n)
{
	int8_t *braid = calloc(min_len + PURE_BRAID_SIZE(1, n), sizeof(int8_t));
	int mark = 0;            // Marker for the braid array index.
	int rand = 0;            // Variable to store RNG values.
	int gen1, gen2;          // Random generators being used.
	int reduced_middle = 0;  // Flag: Remove one reduced middle generator?
	int i;

	// Start with gen2 for better loop management.
	syscall(SYS_getrandom, &rand, sizeof(int), 0);
	gen2 = abs(rand) % (n - 1) + 1;
	for (i = n - 1; i > gen2; i--, mark++)
		braid[mark] = i;
		
	while (true) {
		// Update generators.
		gen1 = gen2;
		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		gen2 = abs(rand) % (n - 1) + 1;

		// Build first generator middle part.
		for (i = 0 + reduced_middle; i < 2; i++, mark++)
			braid[mark] = gen1;

		if (mark >= min_len)
			break;

		// Build freely reduced start and end connecting the two generators.
		if (gen1 >= gen2) {
			for (i = gen1; i > gen2; i--, mark++)
				braid[mark] = i;
			reduced_middle = 0;
		} else if (gen1 < gen2) {
			for (i = gen1 + 1; i < gen2; i++, mark++)
				braid[mark] = i * -1;
			reduced_middle = 1;
		}
	}

	// Build last generator ending.
	for (i = gen1 + 1; i <= n - 1; i++, mark++)
		braid[mark] = i * -1;

	*dest = braid;
	*dest_size = mark;
}

#if 0
/* Generates a pure braid in the following form:
 * high, high-1, ..., low, low, ..., (high-1)^-1, high^-1 */
static void gen_pure_braid(int8_t *dest, int x, int n)
{
	int size = PURE_BRAID_SIZE(x, n);
	int i;

	dest[n - 1 - x] = x;
	dest[n - 1 - x + 1] = x;
	for (i = 0; n > x; i++) {
		dest[i] = n - i;
		dest[size - 1 - i] = -1 * (n - i);
	}
}
#endif

/* Generates a cloaking elemtent in the braid group B_n with the specified delimiters a, b and c, with the specified number of bits in security. */
static void gen_cloaking_elem(int8_t **dest, int *dest_size,
		uint8_t *perm, int a, int b, int n, int sec_level)
{
	// Allocate enough space for the worst case.
	// What about the permutation part at the very start and end?
	int gens_min_len = sec_level / (3 * log2f(n * (n - 1))) + 1;
	int gens_len = 0;
	int8_t *braid = malloc((10 + 4 + 14 * 2 * gens_min_len) * sizeof(int8_t));
	int8_t *pures = NULL;
	int perm_a = 0;
	int perm_b = 0;
	int mark = 0;
	int mark_half;
	int mid;
	int rand;
	int i;

	syscall(SYS_getrandom, &rand, sizeof(int), 0);
	mid = abs(rand) % (n - 1) + 1;

	// Get a and b of inverse permutation.
	for (i = 0; i < n; i++) {
		if (perm[i] == a)
			perm_a = i + 1;
		else if (perm[i] == b)
			perm_b = i + 1;
	}

	// Build permuting part at start to satisfy mid -> a, mid + 1 -> b.
	if (mid <= perm_a && mid < perm_b) {
		for (i = 0; i < perm_a - mid; i++, mark++)
			braid[mark] = perm_a - i - 1;
		if (perm_b < perm_a)
			braid[mark++] = perm_b;
		for (i = 0; i < perm_b - (mid + 1); i++, mark++)
			braid[mark] = perm_b - i - 1;
	} else if (mid < perm_b && mid >= perm_a) {
		for (i = 0; i < perm_b - (mid + 1); i++, mark++)
			braid[mark] = perm_b - i - 1;
		for (i = 0; i < mid - perm_a; i++, mark++)
			braid[mark] = perm_a + i;
	} else if (mid < perm_a && mid >= perm_b) {

		for (i = 0; i < (mid + 1) - perm_b; i++, mark++)
			braid[mark] = perm_b + i;
		if (perm_a != mid + 1) {
			for (i = 0; i < perm_a - mid; i++, mark++)
				braid[mark] = perm_a - i - 1;
			braid[mark++] = mid + 1;
		}
	} else {
		for (i = 0; i < (mid + 1) - perm_b; i++, mark++)
			braid[mark] = perm_b + i;
		if (perm_b < perm_a)
			braid[mark++] = perm_a - 1;
		for (i = 0; i < mid - perm_a; i++, mark++)
			braid[mark] = perm_a + i;
	}

	gen_pure_braids(&pures, &gens_len, gens_min_len, n);
	memcpy(braid + mark, pures, gens_len * sizeof(int8_t));
	mark += gens_len;
	free(pures);

	// Build middle part.
	braid[mark] = mid;
	braid[mark + 1] = mid;
	mark += 2;

	// Build inverse of the first half after the middle.
	mark_half = mark - 2;
	for (i = 0; i < mark_half; i++, mark++)
		braid[mark] = -1 * braid[mark_half - 1 - i];

	*dest = braid;
	*dest_size = mark;
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

#if 0
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

		//gen_pure_braid(braid + mark, xs[i], ys[i], n);
		mark += 2 * (ys[i] - xs[i]) + 2;
	}

	printf("\n\n");

	for (int i = 0; i < size; i++)
		printf("%i, ", braid[i]);
	printf("\n");

	return braid;
}
#endif

