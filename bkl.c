#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "debug.h"

#include "bkl.h"

void left_canonical_form(int8_t **dest, int *dest_len,
	int8_t *braid, int len, int n);

static void inverse(uint8_t *dest, uint8_t *permutation, int n);
static void multiply(uint8_t *dest, uint8_t *perm_a, uint8_t *perm_b, int n);
static bool equal(uint8_t *perm_a, uint8_t *perm_b, int n);
static void perm_to_desc(uint8_t *dest, uint8_t *perm, int n);
static void desc_to_perm(uint8_t *dest, uint8_t *desc, int n);
static void sort_triples(
	uint8_t *dest, uint8_t *a, uint8_t *b, uint8_t *c, int n);
static void meet(uint8_t *dest, uint8_t *desc_a, uint8_t *desc_b, int n);

/* Calculate the inverse of a canonical factor given as a permutation table. */
void inverse(uint8_t *dest, uint8_t *permutation, int n)
{
	uint8_t *temp = malloc(n * sizeof(uint8_t));
	int i;
	for (i = 0; i < n; i++)
		temp[permutation[i] - 1] = i + 1;
	for (i = 0; i < n; i++)
		dest[i] = temp[i];
	free(temp);
}

/* Calculate the product of 2 canonical factors given as permutation tables. */
void multiply(uint8_t *dest, uint8_t *perm_a, uint8_t *perm_b, int n)
{
	uint8_t *temp = malloc(n * sizeof(uint8_t));
	int i;
	for (i = 0; i < n; i++)
		temp[i] = perm_a[perm_b[i] - 1];
	for (i = 0; i < n; i++)
		dest[i] = temp[i];
	free(temp);
}

/* Tests two canonical factors given as permutation tables for equality. */
bool equal(uint8_t *perm_a, uint8_t *perm_b, int n)
{
	int i;
	for (i = 0; i < n; i++) {
		if (perm_a[i] != perm_b[i])
			return false;
	}
	return true;
};

/* Converts a permutation table into a descending cycle decomposition table. */
void perm_to_desc(uint8_t *dest, uint8_t *perm, int n)
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
void desc_to_perm(uint8_t *dest, uint8_t *desc, int n)
{
	uint8_t *temp = malloc(n * sizeof(uint8_t));
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
	free(temp);
}

/* This sorts the triples. I have a feeling that this is kind of a mess and
 * could be improved a lot, especially in terms of simplicity and readability.
 * Additionally, I am not sure if this is the way you are supposed to sort
 * in the first place.
 * If I am going to keep this approach, I should optimize the insertion sort a
 * bit.
 */
void sort_triples(uint8_t *dest, uint8_t *a, uint8_t *b, uint8_t *c, int n)
{
	uint8_t *buckets_a = calloc(n * n, sizeof(uint8_t));
	uint8_t *buckets_b = calloc(n * n, sizeof(uint8_t));
	uint8_t *buckets_c = calloc(n * n, sizeof(uint8_t));
	uint8_t *bucket_i = calloc(n, sizeof(uint8_t));
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
					int temp = *left_b;
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

	free(bucket_i);
	free(buckets_a);
	free(buckets_b);
	free(buckets_c);
}

/* Calculates the meet of two canonical factors A and B given as descending
 * cycle decomposition tables. */
void meet(uint8_t *dest, uint8_t *desc_a, uint8_t *desc_b, int n)
{
	uint8_t *m = malloc(n * sizeof(uint8_t));
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

	free(m);
}

/* Calculates the left-canonical form of a braid given in Artin generators. */
void left_canonical_form(int8_t **dest, int *dest_len,
		int8_t *braid, int len, int n)
{
	int i, j, k;
	int mark = 0;

	uint8_t *trivial = malloc(n * sizeof(uint8_t));
	uint8_t *delta = malloc(n * sizeof(uint8_t));
	uint8_t *delta_inverse = malloc(n * sizeof(uint8_t));
	int8_t *delta_artins = malloc((n - 1) * sizeof(int8_t));
	int8_t *delta_inverse_artins = malloc((n - 1) * sizeof(int8_t));

	uint8_t *band_braid = malloc(len * n * sizeof(uint8_t));
	uint8_t *band_braid_reduced = band_braid;
	uint8_t *band_braid_desc;
	int8_t *artins;

	int *factors_delta_powers = malloc(len * sizeof(int));
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
	free(factors_delta_powers);

	// Convert braid into left canonical form.
	i = 0;
	while (i < len) {
		int t = len;
		for (j = len - 1; j > i; j--) {
			uint8_t *right = band_braid + (j * n);
			uint8_t *left = right - n;
			uint8_t *left_inverse = malloc(n * sizeof(uint8_t));
			uint8_t *left_star = malloc(n * sizeof(uint8_t));
			uint8_t *left_star_desc = malloc(n * sizeof(uint8_t));
			uint8_t *right_desc = malloc(n * sizeof(uint8_t));
			uint8_t *lr_meet_desc = malloc(n * sizeof(uint8_t));
			uint8_t *lr_meet = malloc(n * sizeof(uint8_t));

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

			free(lr_meet);
			free(lr_meet_desc);
			free(right_desc);
			free(left_star_desc);
			free(left_star);
			free(left_inverse);
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
	band_braid_desc = malloc(len * n * sizeof(uint8_t));
	for (i = 0; i < len; i++)
		perm_to_desc(band_braid_desc + i * n, band_braid_reduced + i * n, n);

	free(band_braid);

	// Add fundamental braids as Artin generators to beginning of result.
	// How big can it get?
	artins = malloc(len * n * n * sizeof(int8_t));

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

	free(delta_artins);
	free(delta_inverse_artins);
	
	// Convert rest of braid to Artin generators.
	for (i = 0; i < len; i++) {
		uint8_t *factor = band_braid_desc + (i * n);
		bool *done = calloc(n, sizeof(bool));
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
		free(done);
	}

	free(band_braid_desc);
	free(trivial);
	free(delta);
	free(delta_inverse);

	*dest = artins;
	*dest_len = mark;
}

