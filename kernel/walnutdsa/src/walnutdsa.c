#include <linux/module.h>
#include <linux/types.h>
#include <linux/string.h>
#include <linux/scatterlist.h>
#include <crypto/akcipher.h>
#include <crypto/internal/akcipher.h>

#include "walnut.h"

enum operation {
	SIGN,
	VERIFY
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

union walnut_context {
	enum operation operation;
	struct walnut_sign_data sign_data;
	struct walnut_verify_data verify_data;
};

static int walnut_enc(struct akcipher_request *req);
static int walnut_dec(struct akcipher_request *req);
static int walnut_sign(struct akcipher_request *req);
static int walnut_verify(struct akcipher_request *req);
static int walnut_set_pub_key(struct crypto_akcipher *tfm, const void *data,
	unsigned int data_size);
static int walnut_set_priv_key(struct crypto_akcipher *tfm, const void *data,
	unsigned int data_size);
static unsigned int walnut_max_size(struct crypto_akcipher *tfm);
static void walnut_exit_tfm(struct crypto_akcipher *tfm);

static int free_walnut_context(union walnut_context *context);
static void print_uarr(unsigned int *arr, unsigned int len);
static void print_braid(int8_t *braid, unsigned int len);
static void print_matrix(uint8_t *matrix, unsigned int braid_group);

static void walnut_init(void);
static void walnut_exit(void);


static int walnut_enc(struct akcipher_request *req)
{
	printk(KERN_ALERT "Encryption is not supported by WalnutDSA\n");
	return -1;
}

static int walnut_dec(struct akcipher_request *req)
{
	printk(KERN_ALERT "Decryption is not supported by WalnutDSA\n");
	return -1;
}

static int walnut_sign(struct akcipher_request *req)
{
	struct crypto_akcipher *tfm = crypto_akcipher_reqtfm(req);
	union walnut_context *context = akcipher_tfm_ctx(tfm);

	uint8_t *hash;
	unsigned int hash_size;
	int8_t *sig;
	unsigned int sig_len;

	size_t received;
	int i;

	if (req == NULL) {
		printk(KERN_ALERT "No WalnutDSA request provided\n");
		return -1;
	} else if (context == NULL) {
		printk(KERN_ALERT "No WalnutDSA context provided\n");
		return -1;
	}

	hash_size = req->src_len;
	hash = kmalloc(hash_size * sizeof(uint8_t), 0);
	received = sg_copy_to_buffer(req->src, 1, hash, hash_size);

	if (received == 0) {
		printk(KERN_ALERT "No hash for WalnutDSA signature provided\n");
		return -1;
	}

	generate_signature(&sig, &sig_len, hash, hash_size,
		context->sign_data.braid, context->sign_data.braid_len,
		context->sign_data.a, context->sign_data.b,
		&(context->sign_data.params),
		context->sign_data.rewrite, context->sign_data.sec_level);

	kfree(hash);
	
	sg_copy_from_buffer(req->dst, 1, sig, sig_len * sizeof(int8_t));
	req->dst_len = sig_len * sizeof(int8_t);

	return 0;
}

static int walnut_verify(struct akcipher_request *req)
{
	struct crypto_akcipher *tfm = crypto_akcipher_reqtfm(req);
	union walnut_context *context = akcipher_tfm_ctx(tfm);

	uint8_t *hash;
	unsigned int hash_size;
	int8_t *sig;
	unsigned int sig_len;

	unsigned int received = 0;
	int i;

	if (req == NULL) {
		printk(KERN_ALERT "No WalnutDSA request provided\n");
		return -1;
	} else if (context == NULL) {
		printk(KERN_ALERT "No WalnutDSA context provided\n");
		return -1;
	}

	/*
	 * Data format:
	 * [hash_size][sig_len][hash][sig]
	 *  ...uint..  ..uint.  ....  ...
	 */

	void *buf;
	unsigned int buf_size;

	buf_size = req->src_len;
	buf = kmalloc(buf_size * sizeof(uint8_t), 0);
	received = sg_copy_to_buffer(req->src, 1, buf, buf_size);

	if (received == 0) {
		printk(KERN_ALERT
			"No data for WalnutDSA verification provided\n");
		return -1;
	}

	hash_size = *((unsigned int *) buf);
	sig_len = *((unsigned int *) (buf + sizeof(unsigned int)));

	hash = buf + 2 * sizeof(unsigned int);
	sig = hash + hash_size * sizeof(uint8_t);

	int ret = verify_signature(sig, sig_len, hash, hash_size,
		context->verify_data.matrix, context->verify_data.perm,
		&(context->verify_data.params));

	kfree(buf);
	
	if (ret == -1)
		return -EBADMSG;
	else
		return 0;
}

static int walnut_set_pub_key(struct crypto_akcipher *tfm, const void *data,
	unsigned int data_size)
{
	const struct walnut_verify_data *verify_data = data;
	union walnut_context *context = akcipher_tfm_ctx(tfm);
	void *data_mark;

	if (data == NULL) {
		printk(KERN_ALERT "No verification data provided");
		return -1;
	} else if (free_walnut_context(context)) {
		printk(KERN_ALERT "Could not free previous context");
		return -1;
	}

	context->operation = VERIFY;
	context->verify_data.params = verify_data->params;

	data_mark = data + sizeof(struct walnut_verify_data);
	context->verify_data.params.t_values = kmalloc_array(
		verify_data->params.braid_group, sizeof(unsigned int), 0);
	memcpy(context->verify_data.params.t_values, data_mark,
		verify_data->params.braid_group * sizeof(unsigned int));
	data_mark += verify_data->params.braid_group * sizeof(unsigned int);

	context->verify_data.params.generators
		= kmalloc_array(4, sizeof(unsigned int), 0);
	memcpy(context->verify_data.params.generators, data_mark,
		4 * sizeof(unsigned int));
	data_mark += 4 * sizeof(unsigned int);

	context->verify_data.matrix = kmalloc_array(verify_data->params.braid_group
		* verify_data->params.braid_group, sizeof(uint8_t), 0);
	memcpy(context->verify_data.matrix, data_mark,
		verify_data->params.braid_group
		* verify_data->params.braid_group * sizeof(uint8_t));
	data_mark += verify_data->params.braid_group
		* verify_data->params.braid_group * sizeof(uint8_t);

	context->verify_data.perm = kmalloc_array(
		verify_data->params.braid_group, sizeof(unsigned int), 0);
	memcpy(context->verify_data.perm, data_mark,
		verify_data->params.braid_group * sizeof(unsigned int));

	return 0;
}

static int walnut_set_priv_key(struct crypto_akcipher *tfm, const void *data,
	unsigned int data_size)
{
	struct walnut_sign_data *sign_data = data;
	union walnut_context *context = akcipher_tfm_ctx(tfm);
	void *data_mark;

	if (data == NULL) {
		printk(KERN_ALERT "No signature data provided");
		return -1;
	} else if (free_walnut_context(context)) {
		printk(KERN_ALERT "Could not free previous context");
		return -1;
	}

	context->operation = SIGN;
	context->sign_data.braid_len = sign_data->braid_len;
	context->sign_data.a = sign_data->a;
	context->sign_data.b = sign_data->b;
	context->sign_data.rewrite = sign_data->rewrite;
	context->sign_data.sec_level = sign_data->sec_level;
	context->sign_data.params = sign_data->params;

	data_mark = data + sizeof(struct walnut_sign_data);
	context->sign_data.params.t_values
		= kmalloc_array(sign_data->params.braid_group, sizeof(unsigned int), 0);
	memcpy(context->sign_data.params.t_values, data_mark,
		sign_data->params.braid_group * sizeof(unsigned int));
	data_mark += sign_data->params.braid_group * sizeof(unsigned int);

	context->sign_data.params.generators
		= kmalloc_array(4, sizeof(unsigned int), 0);
	memcpy(context->sign_data.params.generators, data_mark,
		4 * sizeof(unsigned int));
	data_mark += 4 * sizeof(unsigned int);

	context->sign_data.braid
		= kmalloc_array(sign_data->braid_len, sizeof(int8_t), 0);
	memcpy(context->sign_data.braid, data_mark,
		sign_data->braid_len * sizeof(int8_t));

	return 0;
}

static unsigned int walnut_max_size(struct crypto_akcipher *tfm)
{
	/* Giving a maximum size for the signature is not easy.
	 * The current approach is taking all the maximum lengths of the different
	 * parts of the signature and adding them up, which is good.
	 * A problem is the hash size, since it is not known at this point in time.
	 * Thus I chose 256 bytes as a max size for a hash here.
	 * If the BKL normal form is applied, the signature may get much larger.
	 * Dehornoy's reduction could, if it is applied, not reduce much of the
	 * signature, so it is not accounted for here.
	 */

	union walnut_context *context = akcipher_tfm_ctx(tfm);
	if (context == NULL || context->operation == VERIFY) {
		printk(KERN_ALERT
			"No signature data provided for calculating the maximum size. "
			"Run walnut_set_priv_key beforehand.\n");
		return -1;
	}

	unsigned int n = context->sign_data.params.braid_group;
	unsigned int max_len
		= 4 * (10 + 4 + 14 * 2 *
			(context->sign_data.sec_level / (3 * ilog2(n * (n - 1))) + 1))
		+ 2 * context->sign_data.braid_len
		+ 256 * PURE_BRAID_SIZE(1, n) * 2 + 8 + 6;

	if (context->sign_data.rewrite != NONE)
		max_len *= n;

	return max_len * sizeof(int8_t);
}

static void walnut_exit_tfm(struct crypto_akcipher *tfm)
{
	union walnut_context *context = akcipher_tfm_ctx(tfm);
	free_walnut_context(context);
}

static int free_walnut_context(union walnut_context *context)
{
	if (context == NULL)
		return -1;
	
	if (context->operation == SIGN) {
		struct walnut_sign_data data = context->sign_data;
		kfree(data.braid);
		data.braid = NULL;
		kfree(data.params.generators);
		data.params.generators = NULL;
		kfree(data.params.t_values);
		data.params.t_values = NULL;
	} else if (context->operation == VERIFY) {
		struct walnut_verify_data data = context->verify_data;
		kfree(data.perm);
		data.perm = NULL;
		kfree(data.matrix);
		data.matrix = NULL;
		kfree(data.params.generators);
		data.params.generators = NULL;
		kfree(data.params.t_values);
		data.params.t_values = NULL;
	} else {
		return -1;
	}

	return 0;
}

static void print_uarr(unsigned int *arr, unsigned int len) {
	int i;
	char seg[32] = "";
	char *complete = kcalloc(len * 8, sizeof(char), 0);
	for (i = 0; i < len; i++) {
		sprintf(seg, "%u, ", arr[i]);
		strcat(complete, seg);
	}
	printk(KERN_ALERT "%s\n", complete);
	kfree(complete);
}

static void print_braid(int8_t *braid, unsigned int len) {
	int i;
	char seg[5] = "";
	char *complete = kcalloc(len * 5, sizeof(char), 0);
	for (i = 0; i < len; i++) {
		sprintf(seg, "%i, ", braid[i]);
		strcat(complete, seg);
	}
	printk(KERN_ALERT "%s\n", complete);
	kfree(complete);
}

static void print_matrix(uint8_t *matrix, unsigned int braid_group) {
	int i, j;
	char seg[5] = "";
	for (i = 0; i < braid_group; i++) {
		char *complete = kcalloc(braid_group * 5, sizeof(char), 0);
		for (j = 0; j < braid_group; j++) {
			sprintf(seg, "%u, ", matrix[i * braid_group + j]);
			strcat(complete, seg);
		}
		printk(KERN_ALERT "%s\n", complete);
		kfree(complete);
	}
}

static struct akcipher_alg walnut_dsa = {
	.encrypt = walnut_enc,
	.decrypt = walnut_dec,
	.sign = walnut_sign,
	.verify = walnut_verify,
	.set_pub_key = walnut_set_pub_key,
	.set_priv_key = walnut_set_priv_key,
	.max_size = walnut_max_size,
	.exit = walnut_exit_tfm,
	.base = {
		.cra_name = "walnut_dsa",
		.cra_driver_name = "walnut_dsa",
		.cra_priority = 100,
		.cra_module = THIS_MODULE,
		.cra_ctxsize = sizeof(union walnut_context)
	}
};

static int walnut_init(void)
{
	int err = crypto_register_akcipher(&walnut_dsa);

	if (err)
		return err;

	printk(KERN_INFO "WalnutDSA initialized.\n");
	return 0;
}

static void walnut_exit(void)
{
	crypto_unregister_akcipher(&walnut_dsa);
	printk(KERN_INFO "WalnutDSA exited.\n");
}

module_init(walnut_init);
module_exit(walnut_exit);
MODULE_ALIAS_CRYPTO("walnut_dsa");
MODULE_LICENSE("GPL");
MODULE_AUTHOR("Stefan Seil <stefan.seil@campus.lmu.de>");
MODULE_DESCRIPTION("WalnutDSA C implementation");

