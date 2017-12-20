#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <kcapi.h>
#include <errno.h>

#include "galois.h"
#include "debug.h"

enum rewrite_func {
	NONE,
	BKL,
	BKL_DEHORNOY
};

struct pub_params {
	unsigned int braid_group;
	unsigned int galois_order;
	unsigned int *t_values;
	unsigned int *generators;
};

struct walnut_sign_data {
	int8_t *braid;
	unsigned int braid_len;
	unsigned int a;
	unsigned int b;
	unsigned int sec_level;
	enum rewrite_func rewrite;
	struct pub_params params;
};

struct walnut_verify_data {
	uint8_t *matrix;
	unsigned int *perm;
	struct pub_params params;
};

int sign(int8_t **dest, unsigned int *dest_len,
	uint8_t *hash, unsigned int hash_size,
	int8_t *priv_key, unsigned int priv_key_len,
	unsigned int a, unsigned int b,
	struct pub_params *params,
	enum rewrite_func rewrite, unsigned int sec_level);
int verify(int8_t *sig, unsigned int sig_len,
	uint8_t *hash, unsigned int hash_size,
	uint8_t *pub_key_matrix, unsigned int *pub_key_perm,
	struct pub_params *params);
int verify_pub_params(struct pub_params *params);
void gen_priv_key(int8_t **dest, unsigned int *dest_size,
	struct pub_params *params, unsigned int sec_level);
void gen_pub_key(uint8_t **dest_matrix, unsigned int **dest_perm,
	int8_t *priv_key, unsigned int priv_key_len,
	struct pub_params *params);

static unsigned long long fac(unsigned int x);
static void emult(uint8_t *matrix, unsigned int *perm, int artin,
	unsigned int *t_values, unsigned int q, unsigned int n);


int main(void)
{
	int8_t *priv_key;
	unsigned int priv_key_len;
	uint8_t *pub_key_matrix;
	unsigned int *pub_key_perm;
	unsigned int *t_values;
	unsigned int generators[4];
	struct pub_params pub_params;
	unsigned int sec_level;

	unsigned int a, b;
	unsigned int rand;
	int i, j;

	int8_t *signature;
	unsigned int signature_len;

	uint8_t *hash = malloc(32 * sizeof(uint8_t));
	for (i = 0; i < 32; i++) {
		syscall(SYS_getrandom, &rand, sizeof(unsigned int), 0);
		hash[i] = abs(rand) % 256;
	}
	unsigned int hash_size = 32 * sizeof(uint8_t);

	pub_params.braid_group = 8;
	pub_params.galois_order = 256;

	/* Generate random braid_group and galois_order.
	syscall(SYS_getrandom, &rand, sizeof(unsigned int), 0);
	pub_params.braid_group = abs(rand) % (17 - 8) + 8;
	syscall(SYS_getrandom, &rand, sizeof(unsigned int), 0);
	int x = (abs(rand) % (9 - 5) + 5);
	pub_params.galois_order = 1;
	for (i = 0; i < x; i++)
		pub_params.galois_order *= 2;
	*/

	t_values = malloc(pub_params.braid_group * sizeof(unsigned int));
	for (i = 0; i < pub_params.braid_group; i++) {
		syscall(SYS_getrandom, &rand, sizeof(unsigned int), 0);
		t_values[i] = abs(rand) % (pub_params.galois_order - 1) + 1;
	}
	syscall(SYS_getrandom, &rand, sizeof(unsigned int), 0);
	a = abs(rand) % (pub_params.braid_group - 3) + 2;
	t_values[a - 1] = 1;
	if (a == pub_params.braid_group - 2u) {
		b = a + 1;
	} else {
		syscall(SYS_getrandom, &rand, sizeof(unsigned int), 0);
		b = abs(rand) % (pub_params.braid_group - a - 1) + a + 1;
	}
	t_values[b - 1] = 1;
	pub_params.t_values = t_values;

	for (i = 0; i < 4; i++) {
		while (true) {
		again:
			syscall(SYS_getrandom, &rand, sizeof(unsigned int), 0);
			generators[i] = abs(rand) % (pub_params.braid_group - 1) + 1;
			for (j = 0; j < i; j++) {
				if (generators[j] == generators[i])
					goto again;
			}
			break;
		}
	}
	pub_params.generators = generators;

	if (verify_pub_params(&pub_params) != 0) {
		free(t_values);
		return -1;
	}

	sec_level = 128;

	gen_priv_key(&priv_key, &priv_key_len, &pub_params, sec_level);

	gen_pub_key(&pub_key_matrix, &pub_key_perm,
		priv_key, priv_key_len, &pub_params);

	// SIGNATURE GENERATION AND VERIFICATION

	sign(&signature, &signature_len, hash, hash_size, priv_key,
		priv_key_len, a, b, &pub_params, BKL_DEHORNOY, sec_level);

	int ret = verify(signature, signature_len, hash, hash_size,
		pub_key_matrix, pub_key_perm, &pub_params);

	free(pub_key_perm);
	free(pub_key_matrix);
	free(priv_key);
	free(t_values);

	if (ret != 0)
		return -1;
	return 0;
}

int sign(int8_t **dest, unsigned int *dest_len,
	uint8_t *hash, unsigned int hash_size,
	int8_t *priv_key, unsigned int priv_key_len,
	unsigned int a, unsigned int b,
	struct pub_params *params,
	enum rewrite_func rewrite, unsigned int sec_level)
{
	unsigned int data_size = sizeof(struct walnut_sign_data)
		+ sizeof(unsigned int) * params->braid_group
		+ sizeof(unsigned int) * 4
		+ sizeof(int8_t) * priv_key_len;

	void *data = malloc(data_size);
	struct walnut_sign_data *sign_data = data;
	sign_data->braid = NULL;
	sign_data->braid_len = priv_key_len;
	sign_data->a = a;
	sign_data->b = b;
	sign_data->sec_level = sec_level;
	sign_data->rewrite = rewrite;
	sign_data->params = *params;
	sign_data->params.t_values = NULL;
	sign_data->params.generators = NULL;

	void *data_mark = data + sizeof(struct walnut_sign_data);
	memcpy(data_mark, params->t_values,
		params->braid_group * sizeof(unsigned int));
	data_mark += sizeof(unsigned int) * params->braid_group;
	memcpy(data_mark, params->generators, 4 * sizeof(unsigned int));
	data_mark += sizeof(unsigned int) * 4;
	memcpy(data_mark, priv_key, priv_key_len * sizeof(int8_t));

	struct kcapi_handle *handle;
	int err = kcapi_akcipher_init(&handle, "walnut_dsa", 0);
	if (err != 0) {
		char *error_str = strerror(err);
		printf("%s\n", error_str);
	}

	unsigned int max_size = kcapi_akcipher_setkey(
		handle, (uint8_t *) data, data_size);
	if (max_size < 0) {
		perror("kcapi_akcipher_setkey returned with an error");
		return -1;
	}

	int8_t *sig = calloc(max_size, sizeof(int8_t));
	unsigned int sig_len = max_size;

	uint8_t *hash_aligned;
	posix_memalign((void **) &hash_aligned, sysconf(_SC_PAGESIZE), hash_size);
	memcpy(hash_aligned, hash, hash_size);

	unsigned int sig_size = kcapi_akcipher_sign(handle, hash, hash_size,
		(uint8_t *) sig, sig_len * sizeof(int8_t),
		KCAPI_ACCESS_HEURISTIC);
	if (sig_size < 0) {
		perror("kcapi_skcipher_sign return with an error");
	}

	free(hash_aligned);
	kcapi_akcipher_destroy(handle);
	*dest = sig;
	*dest_len = sig_size / sizeof(int8_t);
	return 0;
}

int verify(int8_t *sig, unsigned int sig_len,
		uint8_t *hash, unsigned int hash_size,
		uint8_t *pub_key_matrix, unsigned int *pub_key_perm,
		struct pub_params *params)
{
	unsigned int data_size = sizeof(struct walnut_verify_data)
		+ sizeof(unsigned int) * params->braid_group
		+ sizeof(unsigned int) * 4
		+ sizeof(uint8_t) * params->braid_group * params->braid_group
		+ sizeof(unsigned int) * params->braid_group;

	void *data = malloc(data_size);
	struct walnut_verify_data *sign_data = data;
	sign_data->matrix = NULL;
	sign_data->perm = NULL;
	sign_data->params = *params;
	sign_data->params.t_values = NULL;
	sign_data->params.generators = NULL;

	void *data_mark = data + sizeof(struct walnut_verify_data);
	memcpy(data_mark, params->t_values,
		sizeof(unsigned int) * params->braid_group);
	data_mark += sizeof(unsigned int) * params->braid_group;
	memcpy(data_mark, params->generators, sizeof(unsigned int) * 4);
	data_mark += sizeof(unsigned int) * 4;
	memcpy(data_mark, pub_key_matrix,
		params->braid_group * params->braid_group * sizeof(uint8_t));
	data_mark += sizeof(uint8_t) * params->braid_group * params->braid_group;
	memcpy(data_mark, pub_key_perm, params->braid_group * sizeof(unsigned int));
	data_mark += sizeof(unsigned int) * params->braid_group;

	struct kcapi_handle *handle;
	int err = kcapi_akcipher_init(&handle, "walnut_dsa", 0);
	if (err != 0) {
		char *error_str = strerror(err);
		printf("%s\n", error_str);
	}

	unsigned int max_size = kcapi_akcipher_setpubkey(
		handle, (uint8_t *) data, data_size);
	if (max_size < 0) {
		perror("kcapi_akcipher_setkey returned with an error");
		return -1;
	}

	unsigned int buf_size = 2 * sizeof(unsigned int) + hash_size
		+ sig_len * sizeof(int8_t);
	uint8_t *buf;
	posix_memalign((void **) &buf, sysconf(_SC_PAGESIZE), buf_size);

	unsigned int *buf_mark = (unsigned int *) buf;
	*buf_mark = hash_size;
	buf_mark++;
	*buf_mark = sig_len;
	uint8_t *bfm = buf + 2 * sizeof(unsigned int);
	memcpy(bfm, hash, hash_size);
	bfm += hash_size;
	memcpy(bfm, sig, sig_len * sizeof(int8_t));

	int ret = kcapi_akcipher_verify(handle, buf, buf_size,
		(uint8_t *) sig, max_size, KCAPI_ACCESS_HEURISTIC);
	if (ret == -EBADMSG)
		printf("Verification failed\n");
	else if (ret < 0)
		perror("kcapi_akcipher_verify returned with an error");
	else
		printf("Verification succeeded\n");

	kcapi_akcipher_destroy(handle);

	return ret;
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
void gen_priv_key(int8_t **dest, unsigned int *dest_size,
		struct pub_params *params, unsigned int sec_level)
{
	unsigned int n = params->braid_group;
	float x_prev, x_curr;
	unsigned int len;
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
	len = (unsigned int) x_curr + 1;

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
void gen_pub_key(uint8_t **dest_matrix, unsigned int **dest_perm,
		int8_t *priv_key, unsigned int priv_key_len, struct pub_params *params)
{
	unsigned int n = params->braid_group;
	uint8_t *matrix;
	unsigned int *perm;
	int i;

	matrix = malloc(n * n * sizeof(uint8_t));
	for (i = 0; i < n * n; i++)
		matrix[i] = (i % (n + 1) == 0) ? 1 : 0;

	perm = malloc(n * sizeof(unsigned int));
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
unsigned long long fac(unsigned int x) {
	unsigned long long res = x;
	int i;
	for (i = x - 1; i > 0; i--)
		res *= i;
	return res;
}

static void emult(uint8_t *matrix, unsigned int *perm, int artin,
		unsigned int *t_values, unsigned int q, unsigned int n)
{
	uint8_t *result_matrix = calloc(n * n, sizeof(uint8_t));
	unsigned int *result_perm = calloc(n, sizeof(unsigned int));
	uint8_t artin_matrix[n * n];
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
	for (i = 0; i < n; i++) {
		if (i == abs_artin - 1)
			result_perm[i] = perm[i + 1];
		else if (i == abs_artin + 1 - 1)
			result_perm[i] = perm[i - 1];
		else
			result_perm[i] = perm[i];
	}

	memcpy(matrix, result_matrix, n * n * sizeof(uint8_t));
	memcpy(perm, result_perm, n * sizeof(unsigned int));
	free(result_matrix);
	free(result_perm);
}

