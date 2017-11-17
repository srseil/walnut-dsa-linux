#ifndef WALNUT_H
#define WALNUT_H

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

void generate_signature(int8_t **dest, int *dest_len,
		uint8_t *hash, int hash_size,
		int8_t *priv_key, int priv_key_len,
		uint8_t *pub_key_perm,
		struct pub_params *params, int sec_level);
int verify_signature(int8_t *sig, int sig_len,
		uint8_t *hash, int hash_size,
		uint8_t *pub_key_matrix, uint8_t *pub_key_perm,
		struct pub_params *params);

#endif

