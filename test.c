#define _GNU_SOURCE

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include <stdint.h>

#include <sys/syscall.h>
#include <linux/random.h>

#include "galois.c"
#include "bkl.c"
#include "dehornoy.c"
#include "debug.c"

#define PURE_BRAID_SIZE(x, n) (2 * ((n) - 1 - (x)) + 2)

struct pub_params {
	int braid_group;
	int galois_order;
	uint8_t *t_values;
	int8_t *generators;
	/* Rewriting functions are supposed to be part of the public parameters,
	 * too. Maybe add some flags for which ones to use?
	 */
	/* a and b are supposed to be public parameters, but they are actually only
	 * used during the signature generation, but not during the verification.
	 */
};

/*
void emult(int *matrix, int *perm, int artin, int *t_values, int q, int n);
void gen_cloaking_elem(int **dest, int *dest_size,
	int *perm, int a, int b, int n, int bits);
void gen_pub_key(int **dest_matrix, int **dest_perm,
		int *priv_key, int priv_key_len, int *t_values, int q, int n);
void gen_priv_key(int **dest, int *dest_size, int n, int bits);
void gen_pure_braid(int *dest, int x, int n);
void gen_pure_braids(int **dest, int *dest_size, int min_len, int n);
void encode(int **dest, int *dest_size, char *message, int size, int n);
void generate_signature(int **dest, int *dest_len, char *hash, int hash_size,
		int *priv_key, int priv_key_len, int *pub_key_matrix, int *pub_key_perm,
		int *t_values, int bits);
int verify_signature(int *sig, int sig_len, char *hash, int hash_size,
		int *pub_key_matrix, int *pub_key_perm, int *t_values, int q, int n);
		*/



int verify_pub_params(struct pub_params *params);
void gen_priv_key(int8_t **dest, int *dest_size,
		struct pub_params *params, int bits);
void gen_pub_key(uint8_t **dest_matrix, uint8_t **dest_perm,
		int8_t *priv_key, int priv_key_len, struct pub_params *params);
void generate_signature(int8_t **dest, int *dest_len,
		uint8_t *hash, int hash_size,
		int8_t *priv_key, int priv_key_len,
		uint8_t *pub_key_matrix, uint8_t *pub_key_perm,
		struct pub_params *params, int sec_level);
int verify_signature(int8_t *sig, int sig_len,
		uint8_t *hash, int hash_size,
		uint8_t *pub_key_matrix, uint8_t *pub_key_perm,
		struct pub_params *params);
void encode(int8_t **dest, int *dest_size,
		char *msg, int msg_size,
		int8_t *gens, int n);
void gen_pure_braids(int8_t **dest, int *dest_size, int min_len, int n);
void gen_cloaking_elem(int8_t **dest, int *dest_size,
		uint8_t *perm, int a, int b, int n, int sec_level);
void emult(uint8_t *matrix, uint8_t *perm, int8_t artin,
		uint8_t *t_values, int q, int n);

//void gen_pure_braid(int8_t *dest, int x, int n);

int8_t *clok;
int clok_len;

int main(int argc, char *argv[]) {
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

	// Create public parameters struct.
	struct pub_params params;
	params.braid_group = 8;
	params.galois_order = 32;
	uint8_t t_values[8] = {14, 1, 20, 1, 5, 9, 10, 7};
	//uint8_t t_values[8] = {1, 1, 20, 14, 5, 9, 10, 7};
	params.t_values = t_values;
	int8_t generators[4] = {1, 3, 5, 7};
	params.generators = generators;
	if (verify_pub_params(&params) != 0) {
		return -1;
	}

	int sec_level = 128;

	/*
	int8_t priv_key[] = {
		4, 5, -2, 7, -5, 7, -4, 7, -6, -6, 7, 1, 1, 2, 4, -5, 7, 7, 6, 6, 1, -5,
		7, 5, 3, 1, 5, 2, 5, 7, 5, -4, 5, -2, -4, 2, 4, -1, -1, -6, -4, -1, -2,
		7, -5, 1, -2, 4, 6, -1, -5, 6, 2, -1, 6, -7, -2, -3, -2, -5, -4, -7, -1,
		2, 2, 4, -6, 2, -4, 1, 1, 7, 6, -5, 3, 7, 2, -6, -5, 1, 5, 4, 3, 6, -7,
		-2, 4, 4, -6, 3, -5, 6, 3, 7, -5, -4, -7, -1, -6, -6, 7, -6, -5, -7, 5,
		3, 5, 5, 1, 4, -6, -5, 4, 1
	};
	int priv_key_size = sizeof(priv_key) / sizeof(int8_t);
	*/

	int8_t *priv_key;
	int priv_key_size;
	gen_priv_key(&priv_key, &priv_key_size, &params, sec_level);


	/*
	uint8_t pubk_matrix[] = {
		0, 7, 7, 0, 0, 0, 0, 0,
		0, 7, 7, 14, 14, 0, 0, 0,
		24, 7, 6, 15, 5, 8, 15, 26,
		28, 5, 0, 1, 26, 17, 7, 12,
		6, 11, 0, 30, 19, 17, 6, 13,
		14, 11, 29, 31, 18, 19, 24, 30,
		2, 3, 30, 30, 10, 26, 12, 31,
		0, 0, 0, 0, 0, 0, 0, 1
	};
	uint8_t pubk_perm[] = {3, 6, 1, 5, 2, 7, 8, 4};
	*/

	uint8_t *pubk_matrix, *pubk_perm;
	gen_pub_key(
		&pubk_matrix, &pubk_perm, priv_key, priv_key_size, &params);

	int8_t *signature;
	int signature_len;
	//uint8_t hash[4] = {0x11, 0x22, 0x33, 0x44};
	uint8_t hash[] = {
		0x79, 0xb7, 0xac, 0x30, 0x39, 0x4a, 0xff, 0x82,
		0x92, 0xed, 0x52, 0xf3, 0x8d, 0x52, 0x0c, 0xf8,
		0x2b, 0x01, 0xc8, 0x00, 0x63, 0x07, 0x46, 0xef,
		0x61, 0x4b, 0x26, 0xf0, 0x45, 0xc1, 0x09, 0x5d
	};
	int hash_len = sizeof(hash) / sizeof(uint8_t);


	// Hash "length"? More like size, since it means the bytes?
	generate_signature(&signature, &signature_len, hash, hash_len, priv_key,
		//sizeof(priv_key) / sizeof(int), pubk_matrix, pubk_perm, t_values, 64);
		priv_key_size, pubk_matrix, pubk_perm, &params, sec_level);

	/*
	print_braid(signature, signature_len);
	return 0;
	*/

	verify_signature(signature, signature_len, hash, hash_len,
			pubk_matrix, pubk_perm, &params);

	print_braid(priv_key, priv_key_size);

	puts("");
	print_braid(clok, clok_len);

	return 0;
}

int verify_pub_params(struct pub_params *params) {
	int ones = 0;

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

	for (int i = 0; i < params->braid_group; i++) {
		if (params->t_values[i] < 0
				|| params->t_values[i] >= params->galois_order) {
			// invalid t_value
			return -1;
		}
		ones++;
	}
	if (ones < 2) {
		// t_values are missing ones for a and b
		return -1;
	}

	for (int i = 0; i < 4; i++) {
		if (params->generators[i] < 1
				|| params->generators[i] >= params->braid_group) {
			// invalid generator
			return -1;
		}
		for (int j = 0; j < i; j++) {
			if (params->generators[j] == params->generators[i]) {
				// duplicate generator
				return -1;
			}
		}
	}
	
	return 0;
}

/* Generates a private key in the braid group n with the specified fixed values
 * a and b and the specified security given in bits. */
void gen_priv_key(int8_t **dest, int *dest_size,
		struct pub_params *params, int bits)
{
	int n = params->braid_group;
	int length = 256; // how do we calculate this?
	int8_t *braid = malloc(length * sizeof(int8_t));
	int rand = 0;

	for (int i = 0; i < length; i++) {
		do {
			syscall(SYS_getrandom, &rand, sizeof(int), 0);
			braid[i] = abs(rand) % (2 * n - 1) - (n - 1);
		} while (braid[i] == 0 || i > 0 && braid[i - 1] * -1 == braid[i]);
	}

	*dest = braid;
	*dest_size = length;
}

void gen_pub_key(uint8_t **dest_matrix, uint8_t **dest_perm,
		int8_t *priv_key, int priv_key_len, struct pub_params *params)
{
	int n = params->braid_group;
	uint8_t *matrix = malloc(n * n * sizeof(uint8_t));
	uint8_t *perm = malloc(n * sizeof(uint8_t));

	for (int i = 0; i < n * n; i++)
		matrix[i] = (i % (n + 1) == 0) ? 1 : 0;
	for (int i = 0; i < n; i++)
		perm[i] = i + 1;

	for (int i = 0; i < priv_key_len; i++) {
		emult(matrix, perm, priv_key[i],
			params->t_values, params->galois_order, params->braid_group);
	}

	*dest_matrix = matrix;
	*dest_perm = perm;
}

void generate_signature(int8_t **dest, int *dest_len,
		uint8_t *hash, int hash_size,
		int8_t *priv_key, int priv_key_len,
		uint8_t *pub_key_matrix, uint8_t *pub_key_perm,
		struct pub_params *params, int sec_level)
{
	/* Choose a and b randomly based on the occurencies of 1s in t_values? */
	int n = params->braid_group;

	int a = 0, b = 0;

	int rand = 0; // Helper variable for random number generation.
	int mark = 0; // Helper variable for braid construction.

	int8_t *v, *v1, *v2; // Cloaking elements.
	int v_len, v1_len, v2_len; // Cloaking elements' lengths.

	uint8_t *identity_perm; // Identity permutation for pub key generation.

	int8_t *enc, *sig; // Encoded message and signature.
	int enc_len, sig_len; // Lengths of encoded message and signature.

	identity_perm = malloc(n * sizeof(uint8_t));
	for (int i = 0; i < n; i++)
		identity_perm[i] = i + 1;

	// DOWN HERE!!! i = 0, aber doch i = 1 wegen 1 < a < b < N????
	// Choose a and b based on t_values.
	for (int i = 0; i < n - 1; i++) {
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
	//print_braid(v, v_len);
	//print_braid(v1, v1_len);
	print_braid(v2, v2_len);
	puts("end");

	clok = v2;
	clok_len = v2_len;

	v_len = 0;
	v1_len = 0;
	//v2_len = 0;
	//exit(0);
	// -> v1 und v2 sind schuld

	// Encode hash into braid.
	encode(&enc, &enc_len, hash, hash_size, params->generators, n);
	
	// Build signature.
	sig_len = v_len +  2 * v1_len + v2_len + 2 * priv_key_len + enc_len;
	sig = malloc(sig_len * sizeof(int8_t));
	
	memcpy(sig, v2, v2_len * sizeof(int8_t));
	mark += v2_len;

	for (int i = 0; i < v1_len; i++)
		sig[mark + i] = v1[v1_len - i - 1] * -1;
	mark += v1_len;

	for (int i = 0; i < priv_key_len; i++)
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

	assert(mark == sig_len);

	// Convert signature into left canonical form.
	//left_canonical_form(&sig, &sig_len, sig, sig_len, n);
	printf("SIGNATURE LENGTH: %i (after BKL)\n", sig_len);
	print_braid(sig, sig_len);

	// Apply Dehornoy's algorithm to signature.
	//dehornoy(&sig, &sig_len, sig, sig_len);
	printf("SIGNATURE LENGTH: %i (after Dehornoy)\n", sig_len);
	print_braid(sig, sig_len);

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

	// Encode hash into braid.
	encode(&enc, &enc_len, hash, hash_size, params->generators, n);

	// Generate public key for encoded message.
	gen_pub_key(&msg_pub_key_matrix, &msg_pub_key_perm,
		enc, enc_len, params);


	// Matrix multiplication
	matrix = calloc(n * n, sizeof(uint8_t));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				matrix[i * n + j] ^= galois_mult(msg_pub_key_matrix[i * n + k],
					pub_key_matrix[k * n + j], params->galois_order);
			}
		}
	}

	/*
	print_braid(sig, sig_len);
	exit(0);
	*/

	for (int i = 0; i < sig_len; i++) {
		emult(pub_key_matrix, pub_key_perm, sig[i], params->t_values,
			params->galois_order, params->braid_group);
	}

	printf("Multiplied matrix:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%2i, ", pub_key_matrix[i * n + j]);
		}
		printf("\n");
	}

	printf("Should be equal to:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%2i, ", matrix[i * n + j]);
		}
		printf("\n");
	}
	
	printf("\n");
	for (int i = 0; i < n * n; i++) {
		if (matrix[i] != pub_key_matrix[i]) {
			printf("Verification FAILED! (at %i,%i [%i])\n", i%n, i/n, i);
			return -1;
		}
	}

	printf("Verification SUCCEEDED!\n");
	return 0;
}

/* Encodes a message given as a bit stream into a braid. */
void encode(int8_t **dest, int *dest_size,
		char *msg, int msg_size,
		int8_t *gens, int n)
{
	// Reserve enough space for the unreduced worst case.
	int8_t *braid = calloc(PURE_BRAID_SIZE(1, n) * 4 * msg_size * 2,
		sizeof(int8_t));
	uint8_t digest = *msg;

	int mark = 0;            // Marker for the braid array index.
	//int gens[4] = {0};       // Generators to use in the encoded message.
	int gen1, gen2,          // Generators for the digests, two low bits.
		power1, power2;      // Powers of the generators, two high bits.
	int reduced_middle = 0;  // Flag: Remove one reduced middle generator?
	int digest_complete = 0; // Flag: Get new digest from message?
	int duplicate;           // Flag: Generator already in use?

	/*
	// Pick four random generators to use.
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
		//printf("Randomly chosen generator: %i\n", gens[i]);
	}

	/// WE NEED TO HAVE THOSE GENERATORS FIXED BEFORE!
	gens[0] = 1;
	gens[1] = 3;
	gens[2] = 5;
	gens[3] = 7;
	*/

	// Start with gen2 for better loop management.
	gen2 = gens[(48 & digest) >> 4];
	power2 = ((192 & digest) >> 6) + 1;
	for (int i = n - 1; i > gen2; i--, mark++) {
		braid[mark] = i;
	}
	//printf("Initial: gen: %i, powert: %i\n", gen2, power2);

	// Iterate over the message in 4-bit increments.
	for (int i = 0; i < msg_size * 2 - 1; i++) {
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

		//printf("digest: %u\n", digest);
		//printf("gen: %i, power: %i\n", gen2, power2);

		// Build first generator middle part.
		for (int i = 0 + reduced_middle; i < power1 * 2; i++, mark++) {
			braid[mark] = gen1;
		}

		// Build freely reduced start and end connecting the two generators.
		if (gen1 >= gen2) {
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

void gen_pure_braids(int8_t **dest, int *dest_size, int min_len, int n) {
	int8_t *braid = calloc(min_len + PURE_BRAID_SIZE(1, n), sizeof(int8_t));
	int mark = 0;            // Marker for the braid array index.
	int rand = 0;            // Variable to store RNG values.
	int gen1, gen2;          // Random generators being used.
	int reduced_middle = 0;  // Flag: Remove one reduced middle generator?

	// Start with gen2 for better loop management.
	syscall(SYS_getrandom, &rand, sizeof(int), 0);
	gen2 = abs(rand) % (n - 1) + 1;
	for (int i = n - 1; i > gen2; i--, mark++)
		braid[mark] = i;
		
	while (true) {
		// Update generators.
		gen1 = gen2;
		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		gen2 = abs(rand) % (n - 1) + 1;

		// Build first generator middle part.
		for (int i = 0 + reduced_middle; i < 2; i++, mark++) {
			braid[mark] = gen1;
		}

		if (mark >= min_len)
			break;

		// Build freely reduced start and end connecting the two generators.
		if (gen1 >= gen2) {
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

	// Build last generator ending.
	for (int i = gen1 + 1; i <= n - 1; i++, mark++) {
		braid[mark] = i * -1;
	}

	*dest = braid;
	*dest_size = mark;
}

/* " PRIVATE " */
/* Generates a pure braid in the following form:
 * high, high-1, ..., low, low, ..., (high-1)^-1, high^-1 */
void gen_pure_braid(int8_t *dest, int x, int n) {
	int size = PURE_BRAID_SIZE(x, n);

	dest[n - 1 - x] = x;
	dest[n - 1 - x + 1] = x;
	for (int i = 0; n > x; i++) {
		dest[i] = n - i;
		dest[size - 1 - i] = -1 * (n - i);
	}

	for (int i = 0; i < size; i++)
		printf("%i, ", dest[i]);
	printf("\n");
}

/* Generates a cloaking elemtent in the braid group B_n with the specified delimiters a, b and c, with the specified number of bits in security. */
void gen_cloaking_elem(int8_t **dest, int *dest_size,
		uint8_t *perm, int a, int b, int n, int sec_level) {
	int mark = 0;
	//int gens_len = bits / (2 * log2f((n - 1) * (n - 2))) + 1;
	int gens_len = sec_level / (3 * log2f(n * (n - 1))) + 1;
	printf("Min Generators Length: %i\n", gens_len);

	// Allocate enough space for the worst case.
	// This needs to be refined!
	//int *braid = malloc(PURE_BRAID_SIZE(1, n - 1) * 4 * l * sizeof(int));
	int8_t *braid = malloc((10 + 4 + 14 * 2 * gens_len) * sizeof(int8_t));
	// What about the permutation part at the very start and end???
	
	int rand = 0;
	syscall(SYS_getrandom, &rand, sizeof(int), 0);
	int mid = abs(rand) % (n - 1) + 1;
	printf("Middle choice: %i\n", mid);

	// DEBUG
	printf("DEBUG\n");
	printf("a: %i, b: %i\n", a, b);
	print_braid(perm, n);
	printf("DEBUG\n");

	// Get a and b of inverse permutation.
	int perm_a = 0, perm_b = 0;
	for (int i = 0; i < n; i++) {
		if (perm[i] == a)
			perm_a = i + 1;
		else if (perm[i] == b)
			perm_b = i + 1;
	}
	printf("Perm-A: %i, Perm-B: %i\n", perm_a, perm_b);

	int perm_a_pos = 0, perm_b_pos = 0;
	for (int i = 0; i < n; i++) {
		if (perm[i] == perm_a)
			perm_a_pos = i + 1;
		else if (perm[i] == perm_b)
			perm_b_pos = i + 1;
	}
	printf("Perm-A-Pos: %i, Perm-B-Pos: %i\n", perm_a_pos, perm_b_pos);

	/*
	if (mid <= perm_a_pos) {
		for (int i = 0; i < perm_a_pos - mid; i++, mark++)
			braid[mark] = perm_a_pos - i - 1;
		for (int i = 0; i < perm_b_pos - (mid + 1); i++, mark++)
			braid[mark] = perm_b_pos - i - 1;
	} else if (mid < perm_b_pos) {
		for (int i = 0; i < mid - perm_a_pos; i++, mark++)
			braid[mark] = perm_a_pos + i;
		for (int i = 0; i < perm_b_pos - (mid + 1); i++, mark++)
			braid[mark] = perm_b_pos - i - 1;
	} else {
		for (int i = 0; i < (mid + 1) - perm_b_pos; i++, mark++)
			braid[mark] = perm_b_pos + i;
		for (int i = 0; i < mid - perm_a_pos; i++, mark++)
			braid[mark] = perm_a_pos + i;
	}
	*/

	// Build permuting part at start to satisfy mid -> a, mid + 1 -> b.
	if (mid <= perm_a && mid < perm_b) {
		for (int i = 0; i < perm_a - mid; i++, mark++)
			braid[mark] = perm_a - i - 1;
		if (perm_b < perm_a)
			braid[mark++] = perm_b;
		for (int i = 0; i < perm_b - (mid + 1); i++, mark++)
			braid[mark] = perm_b - i - 1;
	} else if (mid < perm_b && mid >= perm_a) {
		for (int i = 0; i < perm_b - (mid + 1); i++, mark++)
			braid[mark] = perm_b - i - 1;
		for (int i = 0; i < mid - perm_a; i++, mark++)
			braid[mark] = perm_a + i;
	} else if (mid < perm_a && mid >= perm_b) {

		for (int i = 0; i < (mid + 1) - perm_b; i++, mark++)
			braid[mark] = perm_b + i;
		if (perm_a != mid + 1) {
			for (int i = 0; i < perm_a - mid; i++, mark++)
				braid[mark] = perm_a - i - 1;
			braid[mark++] = mid + 1;
		}
		//if (perm_b < perm_a)
		//braid[mark++] = mid + 1;

	} else {
		for (int i = 0; i < (mid + 1) - perm_b; i++, mark++)
			braid[mark] = perm_b + i;
		if (perm_b < perm_a)
			braid[mark++] = perm_a - 1;
		for (int i = 0; i < mid - perm_a; i++, mark++)
			braid[mark] = perm_a + i;
	}

	/*
	// Attach l random pure braid generators.
	for (int i = 0; i < l; i++) {
		int rand = 0, x;
		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		x = abs(rand) % (n - 1) + 1; // 1 <= high < 8
		/*
		do {
			syscall(SYS_getrandom, &rand, sizeof(int), 0);
			low = abs(rand) % (high) + 1; // 1 <= low < high
		} while (low == high);
		printf("Low: %i, High: %i\n", low, high);

		gen_pure_braid(&(braid[mark]), low, high, n);
		mark += PURE_BRAID_SIZE(low, high);
		/

		gen_pure_braid(braid + mark, x, n);
		mark += PURE_BRAID_SIZE(x, n);
	}
	*/

	// Naming is confusing -> size? len? Rathar min and actual size
	int gens_size = 0;
	int8_t *pures = NULL;
	gen_pure_braids(&pures, &gens_size, gens_len, n);
	puts("Generated Pures:");
	print_braid(pures, gens_size);
	memcpy(braid + mark, pures, gens_size * sizeof(int8_t));
	mark += gens_size;

	// Build middle part.
	braid[mark] = mid;
	braid[mark + 1] = mid;
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

void emult(uint8_t *matrix, uint8_t *perm, int8_t artin, uint8_t *t_values, int q, int n) {
	int w = log2(q); // Prime potency for galois field arithmetic
	int abs_artin = abs(artin);
	uint8_t *result_matrix = calloc(n * n, sizeof(uint8_t));
	uint8_t *result_perm = calloc(n, sizeof(uint8_t));

	printf("Artin generator: %i\n", artin);

	// T-variable permutation and evaluation
	uint8_t t_val;
	if (artin < 0) {
		int t_i = perm[(abs_artin - 1 + 1) % n];
		//t_val = galois_inverse(t_values[t_i - 1], w);
		t_val = galois_inverse(t_values[t_i - 1], q);
	} else {
		int t_i = perm[artin - 1];
		t_val = t_values[t_i - 1];
	}

	printf("Evaluated t-value: %i\n", t_val);

	// Colored Burau matrix construction
	uint8_t artin_matrix[n * n];
	for (int row = 0; row < n; row++) {
		if (row == abs_artin - 1) {
			for (int col = 0; col < abs_artin - 2; col++)
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
			for (int col = abs_artin + 1; col < n; col++)
				artin_matrix[row * n + col] = 0;
		} else {
			for (int col = 0; col < n; col++)
				artin_matrix[row * n + col] = (row == col) ? 1 : 0;
		}
	}

	printf("Constructed Colored Burau matrix:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%2i, ", artin_matrix[i * n + j]);
		}
		printf("\n");
	}

	// Matrix multiplication
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				//result_matrix[i * n + j] ^= galois_single_multiply(
				result_matrix[i * n + j] ^= galois_mult(
					matrix[i * n + k], artin_matrix[k * n + j], q);
			}
		}
	}

	printf("Multiplied matrix:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%2i, ", result_matrix[i * n + j]);
		}
		printf("\n");
	}

	// Permutation composition
	for (int i = 0; i < n; i++) {
		if (i == abs_artin - 1)
			result_perm[i] = perm[i + 1];
		else if (i == abs_artin + 1 - 1)
			result_perm[i] = perm[i - 1];
		else
			result_perm[i] = perm[i];
	}

	printf("Composed permutation:\n");
	for (int i = 0; i < n; i++)
		printf("%i ", i + 1);
	printf("\n");
	for (int i = 0; i < n; i++)
		printf("%i ", result_perm[i]);
	printf("\n");

	printf("\n");

	memcpy(matrix, result_matrix, n * n * sizeof(uint8_t));
	memcpy(perm, result_perm, n * sizeof(uint8_t));
	free(result_matrix);
	free(result_perm);
}

/*
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
*/

