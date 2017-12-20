#include <linux/slab.h>
#include <linux/string.h>
#include <linux/kernel.h>

#include "bkl.h"

void left_canonical_form(int8_t **dest, unsigned int *dest_len,
	int8_t *braid, unsigned int len, unsigned int n);

static void inverse(uint8_t *dest, uint8_t *permutation, unsigned int n);
static void multiply(uint8_t *dest, uint8_t *perm_a, uint8_t *perm_b,
		unsigned int n);
static bool equal(uint8_t *perm_a, uint8_t *perm_b, unsigned int n);
static void perm_to_desc(uint8_t *dest, uint8_t *perm, unsigned int n);
static void desc_to_perm(uint8_t *dest, uint8_t *desc, unsigned int n);
static void sort_triples(
	uint8_t *dest, uint8_t *a, uint8_t *b, uint8_t *c, unsigned int n);
static void meet(uint8_t *dest, uint8_t *desc_a, uint8_t *desc_b,
		unsigned int n);

/* Calculate the inverse of a canonical factor given as a permutation table. */
static void inverse(uint8_t *dest, uint8_t *permutation, unsigned int n)
{
	uint8_t *temp = kmalloc(n * sizeof(uint8_t), 0);
	int i;
	for (i = 0; i < n; i++)
		temp[permutation[i] - 1] = i + 1;
	for (i = 0; i < n; i++)
		dest[i] = temp[i];
	kfree(temp);
}

/* Calculate the product of 2 canonical factors given as permutation tables. */
static void multiply(uint8_t *dest, uint8_t *perm_a, uint8_t *perm_b,
		unsigned int n)
{
	uint8_t *temp = kmalloc(n * sizeof(uint8_t), 0);
	int i;
	for (i = 0; i < n; i++)
		temp[i] = perm_a[perm_b[i] - 1];
	for (i = 0; i < n; i++)
		dest[i] = temp[i];
	kfree(temp);
}

/* Tests two canonical factors given as permutation tables for equality. */
static bool equal(uint8_t *perm_a, uint8_t *perm_b, unsigned int n)
{
	int i;
	for (i = 0; i < n; i++) {
		if (perm_a[i] != perm_b[i])
			return false;
	}
	return true;
};

/* Converts a permutation table into a descending cycle decomposition table. */
static void perm_to_desc(uint8_t *dest, uint8_t *perm, unsigned int n)
{
	int i;
	for (i = 0; i < n; i++)
		dest[i] = 0;
	for (i = n - 1; i >= 0; i--) {
		if (dest[i] == 0)
			dest[i] = i + 1;
		if (perm[i] < i + 1)
			dest[perm[i] - 1] = dest[i];
	}
}

/* Converts a descending cycle decomposition table into a permutation table. */
static void desc_to_perm(uint8_t *dest, uint8_t *desc, unsigned int n)
{
	uint8_t *temp = kmalloc(n * sizeof(uint8_t), 0);
	int i;
	for (i = 0; i < n; i++)
		temp[i] = 0;
	for (i = 0; i < n; i++) {
		if (temp[desc[i] - 1] == 0)
			dest[i] = desc[i];
		else
			dest[i] = temp[desc[i] - 1];
		temp[desc[i] - 1] = i + 1;
	}
	kfree(temp);
}

/* This sorts the triples for the meet. This is a mess. This should be rewritten.
 */
static void sort_triples(uint8_t *dest, uint8_t *a, uint8_t *b, uint8_t *c,
		unsigned int n)
{
	uint8_t *buckets_a = kcalloc(n * n, sizeof(uint8_t), 0);
	uint8_t *buckets_b = kcalloc(n * n, sizeof(uint8_t), 0);
	uint8_t *buckets_c = kcalloc(n * n, sizeof(uint8_t), 0);
	uint8_t *bucket_i = kcalloc(n, sizeof(uint8_t), 0);
	int i, j, k;

	for (i = 0; i < n; i++) {
		int bucket = (a[c[i] - 1] - 1) * n;
		int bucket_index = a[c[i] - 1] - 1;
		buckets_a[bucket + bucket_i[bucket_index]] = a[c[i] - 1];
		buckets_b[bucket + bucket_i[bucket_index]] = b[c[i] - 1];
		buckets_c[bucket + bucket_i[bucket_index]] = c[i];
		bucket_i[bucket_index]++;
	}

	for (i = 0; i < n; i++) {
		for (j = 1; j < bucket_i[i]; j++) {
			for (k = j; k > 0; k--) {
				uint8_t *right_b = buckets_b + (i * n + k);
				uint8_t *left_b = buckets_b + (i * n + k - 1);
				uint8_t *right_c = buckets_c + (i * n + k);
				uint8_t *left_c = buckets_c + (i * n + k - 1);

				if (*left_b > *right_b
						|| (*left_b == *right_b && *left_c > *right_c)) {
					uint8_t temp = *left_b;
					*left_b = *right_b;
					*right_b = temp;

					temp = *left_c;
					*left_c = *right_c;
					*right_c = temp;
				} else {
					break;
				}
			}
		}
	}

	j = n - 1;
	j = 0;
	for (i = 0; i < n * n; i++) {
		if (buckets_c[i] != 0)
			dest[j++] = buckets_c[i];
	}

	kfree(bucket_i);
	kfree(buckets_a);
	kfree(buckets_b);
	kfree(buckets_c);
}

/* Calculates the meet of two canonical factors A and B given as descending
 * cycle decomposition tables. */
static void meet(uint8_t *dest, uint8_t *desc_a, uint8_t *desc_b,
		unsigned int n)
{
	uint8_t *m = kmalloc(n * sizeof(uint8_t), 0);
	int i, j;

	for (i = 0; i < n; i++)
		m[i] = n - i;
	sort_triples(m, desc_a, desc_b, m, n);
	j = m[n - 1];
	dest[j - 1] = j;
	for (i = n - 2; i >= 0; i--) {
		if (desc_a[j - 1] != desc_a[m[i] - 1]
				|| desc_b[j - 1] != desc_b[m[i] - 1]) {
			j = m[i];
		}
		dest[m[i] - 1] = j;
	}

	kfree(m);
}

/* Calculates the left-canonical form of a braid given in Artin generators. */
void left_canonical_form(int8_t **dest, unsigned int *dest_len,
		int8_t *braid, unsigned int len, unsigned int n)
{
	int i, j, k;
	unsigned int mark = 0;

	uint8_t *trivial = kmalloc(n * sizeof(uint8_t), 0);
	uint8_t *delta = kmalloc(n * sizeof(uint8_t), 0);
	uint8_t *delta_inverse = kmalloc(n * sizeof(uint8_t), 0);
	int8_t *delta_artins = kmalloc((n - 1) * sizeof(int8_t), 0);
	int8_t *delta_inverse_artins = kmalloc((n - 1) * sizeof(int8_t), 0);

	uint8_t *band_braid = kmalloc(len * n * sizeof(uint8_t), 0);
	uint8_t *band_braid_reduced = band_braid;
	uint8_t *band_braid_desc;
	int8_t *artins;

	int *factors_delta_powers = kmalloc(len * sizeof(int), 0);
	int delta_power = 0;

	for (i = 0; i < n; i++)
		trivial[i] = i + 1;

	delta[0] = n;
	for (i = 1; i < n; i++)
		delta[i] = i;

	delta_inverse[n - 1] = 1;
	for (i = 0; i < n - 1; i++)
		delta_inverse[i] = i + 2;

	for (i = 0; i < len; i++) {
		int abs_braid = abs(braid[i]);
		uint8_t *factor = band_braid + (i * n);

		// Convert into permutation table.
		for (j = 0; j < n; j++) {
			if (j == abs_braid - 1)
				factor[j] = abs_braid + 1;
			else if (j == abs_braid + 1 - 1)
				factor[j] = abs_braid;
			else
				factor[j] = j + 1;
		}

		// Substitue negative generators.
		if (braid[i] < 0) {
			multiply(factor, delta, factor, n);
			factors_delta_powers[i] = -1;
		} else {
			factors_delta_powers[i] = 0;
		}
	}

	for (i = len - 1; i > 0; i--) {
		uint8_t *right = band_braid + (i * n);
		uint8_t *left = right - n;
		int right_delta_power = factors_delta_powers[i];

		if (right_delta_power < 0) {
			// Apply index-shifting automorphism on left factor.
			for (j = 0; j < abs(right_delta_power); j++) {
				multiply(left, delta, left, n);
				multiply(left, left, delta_inverse, n);
			}
			factors_delta_powers[i - 1] += right_delta_power;
		}
	}

	delta_power += factors_delta_powers[0];
	kfree(factors_delta_powers);

	// Convert braid into left canonical form.
	i = 0;
	while (i < len) {
		int t = len;
		for (j = len - 1; j > i; j--) {
			uint8_t *right = band_braid + (j * n);
			uint8_t *left = right - n;
			uint8_t *left_inverse = kmalloc(n * sizeof(uint8_t), 0);
			uint8_t *left_star = kmalloc(n * sizeof(uint8_t), 0);
			uint8_t *left_star_desc = kmalloc(n * sizeof(uint8_t), 0);
			uint8_t *right_desc = kmalloc(n * sizeof(uint8_t), 0);
			uint8_t *lr_meet_desc = kmalloc(n * sizeof(uint8_t), 0);
			uint8_t *lr_meet = kmalloc(n * sizeof(uint8_t), 0);

			inverse(left_inverse, left, n);
			multiply(left_star, left_inverse, delta, n);
			perm_to_desc(left_star_desc, left_star, n);
			perm_to_desc(right_desc, right, n);
			meet(lr_meet_desc, left_star_desc, right_desc, n);
			desc_to_perm(lr_meet, lr_meet_desc, n);

			if (!equal(lr_meet, trivial, n)) {
				t = j - 1;
				multiply(left, left, lr_meet, n);
				inverse(lr_meet, lr_meet, n);
				multiply(right, lr_meet, right, n);
			}

			kfree(lr_meet);
			kfree(lr_meet_desc);
			kfree(right_desc);
			kfree(left_star_desc);
			kfree(left_star);
			kfree(left_inverse);
		}
		i = t + 1;
	}

	// Reduce fundamental braids at the beginning.
	for (; len > 0; len--) {
		if (equal(band_braid_reduced, delta, n)) {
				band_braid_reduced += n;
				delta_power++;
		} else if (equal(band_braid_reduced, delta_inverse, n)) {
				band_braid_reduced += n;
				delta_power--;
		} else {
			break;
		}
	}

	// Reduce trivial braids at the end.
	for (; len > 0
			&& equal(band_braid_reduced + n * (len - 1), trivial, n);
			len--) {
	}

	// Convert band braid into descending cycle decomposition tables.
	band_braid_desc = kmalloc(len * n * sizeof(uint8_t), 0);
	for (i = 0; i < len; i++)
		perm_to_desc(band_braid_desc + i * n, band_braid_reduced + i * n, n);

	kfree(band_braid);

	/* Add fundamental braids as Artin generators to beginning of result.
	 * Is there a better theoretical maximum length for the signature after
	 * the left canonical form? */
	artins = kmalloc(len * n * n * sizeof(int8_t), 0);

	for (i = 0; i < n - 1; i++)
		delta_artins[i] = n - i - 1;
	for (i = 0; i < n - 1; i++)
		delta_inverse_artins[i] = -1 * (i + 1);

	for (i = 0; i < abs(delta_power); i++) {
		if (delta_power > 0)
			memcpy(artins + mark, delta_artins, (n - 1) * sizeof(int8_t));
		else
			memcpy(artins + mark, delta_inverse_artins, (n - 1) * sizeof(int8_t));
		mark += n - 1;
	}

	kfree(delta_artins);
	kfree(delta_inverse_artins);
	
	// Convert rest of braid to Artin generators.
	for (i = 0; i < len; i++) {
		uint8_t *factor = band_braid_desc + (i * n);
		bool *done = kcalloc(n, sizeof(bool), 0);
		for (j = n - 1; j >= 0; j--) {
			int t = j;
			int s = t - 1;

			if (factor[j] !=  j + 1 || done[factor[j] - 1])
				continue;
			
			for (; s >= 0; s--) {
				if (factor[s] != factor[t])
					continue;
				for (k = 0; k < t - s; k++)
					artins[mark++] = t - k;
				for (k = 0; k < t - s - 1; k++)
					artins[mark++] = -1 * (s + k + 2);
				t = s;
			}
			done[factor[t] - 1] = true;
		}
		kfree(done);
	}

	kfree(band_braid_desc);
	kfree(trivial);
	kfree(delta);
	kfree(delta_inverse);

	*dest = artins;
	*dest_len = mark;
}

