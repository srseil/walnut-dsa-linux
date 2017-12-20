#include <linux/module.h>
#include <linux/slab.h>
#include <linux/types.h>
#include <linux/random.h>
#include <linux/string.h>
#include <linux/scatterlist.h>
#include <crypto/akcipher.h>
#include <linux/log2.h>
#include <linux/kernel.h>

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
void gen_priv_key(int8_t **dest, unsigned int *dest_size,
	struct pub_params *params, unsigned int sec_level);
void gen_pub_key(uint8_t **dest_matrix, unsigned int **dest_perm,
	int8_t *priv_key, unsigned int priv_key_len,
	struct pub_params *params);

static unsigned long long fac(unsigned int x);
static void emult(uint8_t *matrix, unsigned int *perm, int artin,
	unsigned int *t_values, unsigned int q, unsigned int n);
static uint8_t galois_mult(unsigned int a, unsigned int b, unsigned int q);
static uint8_t galois_inverse(unsigned int x, unsigned int q);

static void print_braid(int8_t *braid, unsigned int len);
static int test_init(void);
static void test_exit(void);


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

	uint8_t *hash = kmalloc(32 * sizeof(uint8_t), 0);
	for (i = 0; i < 32; i++) {
		get_random_bytes(&rand, sizeof(unsigned int));
		hash[i] = abs(rand) % 256;
	}
	unsigned int hash_size = 32 * sizeof(uint8_t);

	pub_params.braid_group = 8;
	pub_params.galois_order = 256;

	/* Generate random braid_group and galois_order.
	get_random_bytes(&rand, sizeof(unsigned int));
	pub_params.braid_group = abs(rand) % (17 - 8) + 8;
	get_random_bytes(&rand, sizeof(unsigned int));
	int x = (abs(rand) % (9 - 5) + 5);
	pub_params.galois_order = 1;
	for (i = 0; i < x; i++)
		pub_params.galois_order *= 2;
	*/

	t_values = kmalloc(pub_params.braid_group * sizeof(unsigned int), 0);
	for (i = 0; i < pub_params.braid_group; i++) {
		get_random_bytes(&rand, sizeof(unsigned int));
		t_values[i] = abs(rand) % (pub_params.galois_order - 1) + 1;
	}
	get_random_bytes(&rand, sizeof(unsigned int));
	a = abs(rand) % (pub_params.braid_group - 3) + 2;
	t_values[a - 1] = 1;
	if (a == pub_params.braid_group - 2u) {
		b = a + 1;
	} else {
		get_random_bytes(&rand, sizeof(unsigned int));
		b = abs(rand) % (pub_params.braid_group - a - 1) + a + 1;
	}
	t_values[b - 1] = 1;
	pub_params.t_values = t_values;

	for (i = 0; i < 4; i++) {
		while (true) {
		again:
			get_random_bytes(&rand, sizeof(unsigned int));
			generators[i] = abs(rand) % (pub_params.braid_group - 1) + 1;
			for (j = 0; j < i; j++) {
				if (generators[j] == generators[i])
					goto again;
			}
			break;
		}
	}
	pub_params.generators = generators;

	sec_level = 128;

	gen_priv_key(&priv_key, &priv_key_len, &pub_params, sec_level);

	gen_pub_key(&pub_key_matrix, &pub_key_perm,
		priv_key, priv_key_len, &pub_params);

	// SIGNATURE GENERATION AND VERIFICATION
	sign(&signature, &signature_len, hash, hash_size, priv_key,
		priv_key_len, a, b, &pub_params, BKL_DEHORNOY, sec_level);

	int ret = verify(signature, signature_len, hash, hash_size,
		pub_key_matrix, pub_key_perm, &pub_params);

	kfree(pub_key_perm);
	kfree(pub_key_matrix);
	kfree(priv_key);
	kfree(t_values);

	return ret;
}

int verify(int8_t *sig, unsigned int sig_len,
		uint8_t *hash, unsigned int hash_size,
		uint8_t *pub_key_matrix, unsigned int *pub_key_perm,
		struct pub_params *params)
{
	int err;
	unsigned int max_size;
	struct crypto_akcipher *tfm;

	unsigned int data_size = sizeof(struct walnut_verify_data)
		+ sizeof(unsigned int) * params->braid_group
		+ sizeof(unsigned int) * 4
		+ sizeof(uint8_t) * params->braid_group * params->braid_group
		+ sizeof(unsigned int) * params->braid_group;
	uint8_t *data = kmalloc(data_size, 0);
	struct walnut_verify_data *sign_data = data;
	sign_data->matrix = NULL;
	sign_data->perm = NULL;
	sign_data->params = *params;
	sign_data->params.t_values = NULL;
	sign_data->params.generators = NULL;

	uint8_t *data_mark = data + sizeof(struct walnut_verify_data);
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

	tfm = crypto_alloc_akcipher("walnut_dsa", CRYPTO_ALG_TYPE_AKCIPHER, 0);
	if (IS_ERR(tfm)) {
		printk(KERN_ALERT "ERROR: crypto_alloc_akcipher\n");
		return;
	}

	err = crypto_akcipher_set_pub_key(tfm, data, data_size);
	if (err != 0)
		printk(KERN_ALERT "ERROR: crypto_akcipher_set_pub_key\n");

	kfree(data);

	struct scatterlist sg_src[1], sg_dst[1];
	sg_init_table(sg_src, 1);
	sg_mark_end(sg_src);
	sg_init_table(sg_dst, 1);
	sg_mark_end(sg_dst);

	unsigned int buf_size = 2 * sizeof(unsigned int) + hash_size
		+ sig_len * sizeof(int8_t);
	uint8_t *buf = kmalloc(buf_size, 0);

	*((unsigned int *) buf) = hash_size;
	*((unsigned int *) (buf + sizeof(unsigned int))) = sig_len;
	memcpy(buf + 2 * sizeof(unsigned int), hash, hash_size);
	memcpy(buf + 2 * sizeof(unsigned int) + hash_size,
		sig, sig_len * sizeof(int8_t));
	sg_set_buf(sg_src, buf, buf_size);

	struct akcipher_request *req = akcipher_request_alloc(tfm, 0);
	akcipher_request_set_crypt(req, sg_src, sg_dst, buf_size, 0);

	int ret = crypto_akcipher_verify(req);

	kfree(buf);

	if (ret == -EBADMSG)
		printk(KERN_ALERT "Verification failed");
	else if (ret < 0)
		printk(KERN_ALERT "Error code: %i", ret);
	else
		printk(KERN_ALERT "Verification succeeded");

	crypto_free_akcipher(tfm);

	return ret;
}

int sign(int8_t **dest, unsigned int *dest_len,
	uint8_t *hash, unsigned int hash_size,
	int8_t *priv_key, unsigned int priv_key_len,
	unsigned int a, unsigned int b,
	struct pub_params *params,
	enum rewrite_func rewrite, unsigned int sec_level)
{
	int err;
	unsigned int max_size;
	struct crypto_akcipher *tfm;
	unsigned int data_size = sizeof(struct walnut_sign_data)
		+ sizeof(unsigned int) * params->braid_group
		+ sizeof(unsigned int) * 4
		+ sizeof(int8_t) * priv_key_len;

	uint8_t *data = kmalloc(data_size, 0);
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

	uint8_t *data_mark = data + sizeof(struct walnut_sign_data);
	memcpy(data_mark, params->t_values,
		params->braid_group * sizeof(unsigned int));
	data_mark += sizeof(unsigned int) * params->braid_group;
	memcpy(data_mark, params->generators, 4 * sizeof(unsigned int));
	data_mark += sizeof(unsigned int) * 4;
	memcpy(data_mark, priv_key, priv_key_len * sizeof(int8_t));
	

	tfm = crypto_alloc_akcipher("walnut_dsa", CRYPTO_ALG_TYPE_AKCIPHER, 0);
	if (IS_ERR(tfm)) {
		printk(KERN_ALERT "ERROR: crypto_alloc_akcipher\n");
		return -1;
	}

	err = crypto_akcipher_set_priv_key(tfm, data, data_size);
	if (err != 0)
		printk(KERN_ALERT "ERROR: crypto_akcipher_set_priv_key\n");
	kfree(data);

	max_size = crypto_akcipher_maxsize(tfm);
	
	struct scatterlist sg_src[1], sg_dst[1];
	sg_init_table(sg_src, 1);
	sg_mark_end(sg_src);
	sg_init_table(sg_dst, 1);
	sg_mark_end(sg_dst);

	sg_set_buf(sg_src, hash, hash_size);

	int8_t *sig = kcalloc(max_size, sizeof(int8_t), 0);
	unsigned int sig_len = max_size * sizeof(int8_t);
	sg_set_buf(sg_dst, sig, sig_len);

	int i;

	struct akcipher_request *req = akcipher_request_alloc(tfm, 0);
	akcipher_request_set_crypt(req, sg_src, sg_dst, hash_size, sig_len);

	crypto_akcipher_sign(req);

	*dest = sig;
	*dest_len = req->dst_len;

	akcipher_request_free(req);
	crypto_free_akcipher(tfm);
	
	return 0;
}

static void print_braid(int8_t *braid, unsigned int len) {
	int i;
	char *seg = kcalloc(8, sizeof(char), 0);
	char *complete = kcalloc(len * 8, sizeof(char), 0);

	for (i = 0; i < len; i++) {
		sprintf(seg, "%i, ", braid[i]);
		strncat(complete, seg, 8 * sizeof(char));
	}

	printk(KERN_ALERT "%s\n", complete);
	kfree(seg);
	kfree(complete);
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

	/* A GENERAL LOG2 IS NEEDED FOR THIS, BUT IT IS NOT AVAILABLE?
	// Get minimum key length with Newton's method.
	x_prev = sec_level;
	for (i = 0; i < 3; i++) {
		x_curr = x_prev
			- (x_prev + (n - 2) * ilog2(x_prev) - (sec_level + ilog2(fac(n - 1))))
			/ ((n - 2) / (x_prev * ilog(2)) + 1);
		x_prev = x_curr;
	}
	len = (unsigned int) x_curr + 1;
	*/

	if (sec_level == 112)
		len = 86;
	else if (sec_level == 128)
		len = 101;
	else if (sec_level == 192)
		len = 161;
	else if (sec_level == 256)
		len = 222;
	else
		len = sec_level;

	braid = kmalloc(len * sizeof(int8_t), 0);
	for (i = 0; i < len; i++) {
		do {
			get_random_bytes(&rand, sizeof(int));
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
	uint8_t *result_matrix = kcalloc(n * n, sizeof(uint8_t), 0);
	unsigned int *result_perm = kcalloc(n, sizeof(unsigned int), 0);
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
	kfree(result_matrix);
	kfree(result_perm);
}

/* Irreducable polynomials for the finite fields GF(2^n) with n = 5, 6, 7, 8. */
static unsigned int irr_polys[4] = {
	37,  /* x^5 + x^2 + 1 */
	67,  /* x^6 + x + 1 */
	131, /* x^7 + x + 1 */
	283, /* x^8 + x^4 + x^3 + x + 1 */
};

static uint8_t galois_mult(unsigned int a, unsigned int b, unsigned int q) {
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


static uint8_t galois_inverse(unsigned int x, unsigned int q) {
	int i;
	for (i = 0; i < q; i++) {
		if (galois_mult(x, i, q) == 1)
			return i;
	}
	return 0;
}

static int test_init(void)
{
	printk(KERN_ALERT "Walnut Test starting...\n");

	main();

	return 0;
}

static void test_exit(void)
{
	printk(KERN_ALERT "Walnut Test exiting...\n");
}

module_init(test_init);
module_exit(test_exit);
MODULE_LICENSE("GPL");

