#ifndef WALNUT_H
#define WALNUT_H

#define PURE_BRAID_SIZE(x, n) (2 * ((n) - 1 - (x)) + 2)

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

#endif

