#include <linux/slab.h>
#include <linux/types.h>
#include <linux/string.h>
#include <linux/kernel.h>
#include <linux/log2.h>
#include <linux/random.h>

#include "galois.h"
#include "bkl.h"
#include "dehornoy.h"

#include "walnut.h"

void generate_signature(int8_t **dest, unsigned int *dest_len,
		uint8_t *hash, unsigned int hash_size,
		int8_t *priv_key, unsigned int priv_key_len,
		unsigned int a, unsigned int b,
		struct pub_params *params,
		enum rewrite_func rewrite, unsigned int sec_level);
int verify_signature(int8_t *sig, unsigned int sig_len,
		uint8_t *hash, unsigned int hash_size,
		uint8_t *pub_key_matrix, unsigned int *pub_key_perm,
		struct pub_params *params);

static void gen_pub_key(uint8_t **dest_matrix, unsigned int **dest_perm,
		int8_t *priv_key, unsigned int priv_key_len, struct pub_params *params);
static void encode(int8_t **dest, unsigned int *dest_size,
		uint8_t *msg, unsigned int msg_size,
		unsigned int *gens, unsigned int n);
static void gen_pure_braids(int8_t **dest, unsigned int *dest_size,
		unsigned int min_len, unsigned int n);
static void gen_cloaking_elem(int8_t **dest, unsigned int *dest_size,
		unsigned int *perm, unsigned int a, unsigned int b,
		unsigned int n, unsigned int sec_level);
static void emult(uint8_t *matrix, unsigned int *perm, int artin,
		unsigned int *t_values, unsigned int q, unsigned int n);


static void gen_pub_key(uint8_t **dest_matrix, unsigned int **dest_perm,
		int8_t *priv_key, unsigned int priv_key_len, struct pub_params *params)
{
	unsigned int n = params->braid_group;
	uint8_t *matrix;
	unsigned int *perm;
	int i;

	matrix = kmalloc(n * n * sizeof(uint8_t), 0);
	for (i = 0; i < n * n; i++)
		matrix[i] = (i % (n + 1) == 0) ? 1 : 0;

	perm = kmalloc(n * sizeof(unsigned int), 0);
	for (i = 0; i < n; i++)
		perm[i] = i + 1;

	for (i = 0; i < priv_key_len; i++) {
		emult(matrix, perm, priv_key[i],
			params->t_values, params->galois_order, params->braid_group);
	}

	*dest_matrix = matrix;
	*dest_perm = perm;
}

void generate_signature(int8_t **dest, unsigned int *dest_len,
		uint8_t *hash, unsigned int hash_size,
		int8_t *priv_key, unsigned int priv_key_len,
		unsigned int a, unsigned int b,
		struct pub_params *params,
		enum rewrite_func rewrite, unsigned int sec_level)
{
	int8_t *v, *v1, *v2; // Cloaking elements.
	unsigned int v_len, v1_len, v2_len; // Cloaking elements' lengths.

	unsigned int *identity_perm; // Identity permutation for cloaking element.
	unsigned int *key_perm; // Key permutation for cloaking elements.

	int8_t *enc, *sig, *sig_old; // Encoded message and signature.
	unsigned int enc_len, sig_len; // Lengths of encoded message and signature.

	unsigned int n = params->braid_group;
	unsigned int mark = 0;
	int i;

	identity_perm = kmalloc(n * sizeof(unsigned int), 0);
	for (i = 0; i < n; i++)
		identity_perm[i] = i + 1;

	key_perm = kmalloc(n * sizeof(unsigned int), 0);
	for (i = 0; i < n; i++)
		key_perm[i] = i + 1;
	for (i = 0; i < priv_key_len; i++) {
		int abs_artin = abs(priv_key[i]);
		unsigned int tmp = key_perm[abs_artin - 1];
		key_perm[abs_artin - 1] = key_perm[abs_artin];
		key_perm[abs_artin] = tmp;
	}

	// Generate cloaking elements.
	gen_cloaking_elem(&v, &v_len, identity_perm, a, b, n, sec_level);
	gen_cloaking_elem(&v1, &v1_len, key_perm, a, b, n, sec_level);
	gen_cloaking_elem(&v2, &v2_len, key_perm, a, b, n, sec_level);
	kfree(key_perm);
	kfree(identity_perm);

	// Encode hash into braid.
	encode(&enc, &enc_len, hash, hash_size, params->generators, n);
	
	// Build signature.
	sig_len = v_len +  2 * v1_len + v2_len + 2 * priv_key_len + enc_len;
	sig = kmalloc(sig_len * sizeof(int8_t), 0);
	
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

	kfree(enc);
	kfree(v2);
	kfree(v1);
	kfree(v);

	// Convert signature into left canonical form.
	if (rewrite == BKL || rewrite == BKL_DEHORNOY) {
		sig_old = sig;
		left_canonical_form(&sig, &sig_len, sig, sig_len, n);
		kfree(sig_old);
	}

	// Apply Dehornoy's algorithm to signature.
	if (rewrite == BKL_DEHORNOY)
		dehornoy(&sig, &sig_len, sig, sig_len);

	*dest = sig;
	*dest_len = sig_len;
}

int verify_signature(int8_t *sig, unsigned int sig_len,
		uint8_t *hash, unsigned int hash_size,
		uint8_t *pub_key_matrix, unsigned int *pub_key_perm,
		struct pub_params *params)
{
	unsigned int n = params->braid_group;
	int8_t *enc;
	unsigned int enc_len;
	uint8_t *msg_pub_key_matrix;
	unsigned int *msg_pub_key_perm;
	uint8_t *matrix;
	int i, j, k;

	// Encode hash into braid.
	encode(&enc, &enc_len, hash, hash_size, params->generators, n);

	// Generate public key for encoded message.
	gen_pub_key(&msg_pub_key_matrix, &msg_pub_key_perm,
		enc, enc_len, params);
	kfree(enc);

	// Matrix multiplication
	matrix = kcalloc(n * n, sizeof(uint8_t), 0);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				matrix[i * n + j] ^= galois_mult(msg_pub_key_matrix[i * n + k],
					pub_key_matrix[k * n + j], params->galois_order);
			}
		}
	}
	kfree(msg_pub_key_matrix);
	kfree(msg_pub_key_perm);

	for (i = 0; i < sig_len; i++) {
		emult(pub_key_matrix, pub_key_perm, sig[i], params->t_values,
			params->galois_order, params->braid_group);
	}

	for (i = 0; i < n * n; i++) {
		if (matrix[i] != pub_key_matrix[i]) {
			kfree(matrix);
			return -1;
		}
	}

	kfree(matrix);
	return 0;
}

/* Encodes a message given as a bit stream into a braid. */
static void encode(int8_t **dest, unsigned int *dest_size,
		uint8_t *msg, unsigned int msg_size,
		unsigned int *gens, unsigned int n)
{
	// Reserve enough space for the unreduced worst case.
	int8_t *braid = kcalloc(
		msg_size * PURE_BRAID_SIZE(1, n) * 2 + 8 + 6, sizeof(int8_t), 0);
	uint8_t digest = msg[0]; // Currently encoded hash digest.
	unsigned int mark = 0;   // Marker for the braid array index.
	unsigned int gen1, gen2, // Generators for the digests, two low bits.
		power1, power2;      // Powers of the generators, two high bits.
	int reduced_middle = 0;  // Flag: Remove one reduced middle generator?
	int digest_complete = 0; // Flag: Get new digest from message?
	int i, j;

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

static void gen_pure_braids(int8_t **dest, unsigned int *dest_size,
		unsigned int min_len, unsigned int n)
{
	int8_t *braid = kcalloc(min_len + PURE_BRAID_SIZE(1, n), sizeof(int8_t), 0);
	unsigned int rand = 0;   // Variable to store RNG values.
	unsigned int gen1, gen2; // Random generators being used.
	int reduced_middle = 0;  // Flag: Remove one reduced middle generator?
	unsigned int mark = 0;
	int i;

	// Start with gen2 for better loop management.
	get_random_bytes(&rand, sizeof(unsigned int));
	gen2 = abs(rand) % (n - 1) + 1;
	for (i = n - 1; i > gen2; i--, mark++)
		braid[mark] = i;
		
	while (true) {
		// Update generators.
		gen1 = gen2;
		get_random_bytes(&rand, sizeof(unsigned int));
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

/* Generates a cloaking elemtent in the braid group B_n with the specified delimiters a, b and c, with the specified number of bits in security. */
static void gen_cloaking_elem(int8_t **dest, unsigned int *dest_size,
		unsigned int *perm, unsigned int a, unsigned int b,
		unsigned int n, unsigned int sec_level)
{
	unsigned int gens_min_len = sec_level / (3 * ilog2(n * (n - 1))) + 1;
	unsigned int gens_len = 0;
	int8_t *braid = kcalloc(10 + 4 + 14 * 2 * gens_min_len, sizeof(int8_t), 0);
	int8_t *pures = NULL;
	unsigned int perm_a = 0;
	unsigned int perm_b = 0;
	unsigned int mark = 0;
	unsigned int mark_half;
	unsigned int mid;
	unsigned int rand;
	int i;

	get_random_bytes(&rand, sizeof(unsigned int));
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
	kfree(pures);

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

static void emult(uint8_t *matrix, unsigned int *perm, int artin,
		unsigned int *t_values, unsigned int q, unsigned int n)
{
	uint8_t *result_matrix = kcalloc(n * n, sizeof(uint8_t), 0);
	uint8_t *artin_matrix = kcalloc(n * n, sizeof(uint8_t), 0);
	unsigned int t_val;
	int abs_artin = abs(artin);
	int row, col;
	int i, j, k;

	// T-variable permutation and evaluation
	if (artin < 0) {
		unsigned int t_i = perm[(abs_artin - 1 + 1) % n];
		t_val = galois_inverse(t_values[t_i - 1], q);
	} else {
		unsigned int t_i = perm[artin - 1];
		t_val = t_values[t_i - 1];
	}

	// Colored Burau matrix construction
	for (row = 0; row < n; row++) {
		if (row == abs_artin - 1) {
			for (col = 0; col < abs_artin - 2; col++)
				artin_matrix[row * n + col] = 0;
			if (artin < 0) {
				if (artin != -1)
					artin_matrix[row * n + abs_artin - 2] = 1;
				artin_matrix[row * n + abs_artin - 1] = t_val;
				if (abs_artin != n)
					artin_matrix[row * n + abs_artin] = t_val;
			} else {
				if (artin != 1)
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

	// Matrix multiplication
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				result_matrix[i * n + j] ^= galois_mult(
					matrix[i * n + k], artin_matrix[k * n + j], q);
			}
		}
	}

	// Permutation composition
	unsigned int tmp = perm[abs_artin - 1];
	perm[abs_artin - 1] = perm[abs_artin];
	perm[abs_artin] = tmp;

	memcpy(matrix, result_matrix, n * n * sizeof(uint8_t));
	kfree(artin_matrix);
	kfree(result_matrix);
}

