#ifndef BKL
#define BKL

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define NONE 0
#define POSITIVE 1
#define NEGATIVE 2

void left_canonical_form(int8_t **dest, int *dest_len, int8_t *braid, int len, int n);

/*
struct band_gen {
	int t;
	int s;
};
*/

/*
struct canonical_factor {
	int fundamental;       // Flag for fundamental braid (none, pos, neg).
	int fundamental_power; // The power of the fundamental braid.
	//int *generators;       // Pointer to array of Artin generators.
	struct band_gen *generators;
	int num_generators;    // Number of Artin generators in this factor.
};
*/

/* Calculate the inverse of a canonical factor given in a permutation table. */
void inverse(uint8_t *dest, uint8_t *permutation, int n) {
	uint8_t *temp = malloc(n * sizeof(uint8_t));
	for (int i = 0; i < n; i++)
		temp[permutation[i] - 1] = i + 1;
	for (int i = 0; i < n; i++)
		dest[i] = temp[i];
	free(temp);
}

/* Calculate the product of 2 canonical factors given in permutation tables. */
void multiply(uint8_t *dest, uint8_t *perm_a, uint8_t *perm_b, int n) {
	uint8_t *temp = malloc(n * sizeof(uint8_t));
	for (int i = 0; i < n; i++)
		temp[i] = perm_a[perm_b[i] - 1];
	for (int i = 0; i < n; i++)
		dest[i] = temp[i];
	free(temp);
}

/* Tests two canonical factors given as permutation tables for equality. */
bool equal(uint8_t *perm_a, uint8_t *perm_b, int n) {
	for (int i = 0; i < n; i++) {
		if (perm_a[i] != perm_b[i])
			return false;
	}
	return true;
};

/* Converts a permutation table into a descending cycle decomposition table. */
void perm_to_desc(uint8_t *dest, uint8_t *perm, int n) {
	for (int i = 0; i < n; i++)
		dest[i] = 0;
	for (int i = n - 1; i >= 0; i--) {
		if (dest[i] == 0)
			dest[i] = i + 1;
		if (perm[i] < i + 1)
			dest[perm[i] - 1] = dest[i];
	}
}

/* Converts a descending cycle decomposition table into a permutation table. */
void desc_to_perm(uint8_t *dest, uint8_t *desc, int n) {
	uint8_t *temp = malloc(n * sizeof(uint8_t));
	for (int i = 0; i < n; i++)
		temp[i] = 0;
	for (int i = 0; i < n; i++) {
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
void sort_triples(uint8_t *dest, uint8_t *a, uint8_t *b, uint8_t *c, int n) {
	// Use bucket sort with insertion sort.
	uint8_t *buckets_a = calloc(n * n, sizeof(uint8_t));
	uint8_t *buckets_b = calloc(n * n, sizeof(uint8_t));
	uint8_t *buckets_c = calloc(n * n, sizeof(uint8_t));


	/*
	int a[8] = {5, 3, 3, 5, 5, 6, 7, 8};
	int b[8] = {4, 4, 3, 4, 5, 6, 7, 8};
	// cs[8] = {8, 7, 6, 5, 4, 3, 2, 1};
	//
	//          333, 342, 541, 544, 555, 666, 777, 888
	meet(result, a, b, 8);
	*/

	/*
	printf("\tTEST:\n");
	for (int i = 0; i < n; i++) {
		printf("%i ", b[i]);
	}
	printf("\n");
	*/


	uint8_t *bucket_i = calloc(n, sizeof(uint8_t));
	for (int i = 0; i < n; i++) {
		int bucket = (a[c[i] - 1] - 1) * n;
		int bucket_index = a[c[i] - 1] - 1;
		buckets_a[bucket + bucket_i[bucket_index]] = a[c[i] - 1];
		buckets_b[bucket + bucket_i[bucket_index]] = b[c[i] - 1];
		buckets_c[bucket + bucket_i[bucket_index]] = c[i];
		bucket_i[bucket_index]++;
		//printf("%i: %i\n", bucket_index, bucket_i[bucket_index]);
	}

	for (int i = 0; i < n; i++) {
		for (int j = 1; j < bucket_i[i]; j++) {
			//printf("Max: %i\n", bucket_i[i]);
			for (int k = j; k > 0; k--) {
				uint8_t *right_b = buckets_b + (i * n + k);
				uint8_t *left_b = buckets_b + (i * n + k - 1);
				uint8_t *right_c = buckets_c + (i * n + k);
				uint8_t *left_c = buckets_c + (i * n + k - 1);

				//printf("Comparison: %i %i ... %i %i", *left_b, *right_b, *left_c, *right_c);

				if (*left_b > *right_b
						|| *left_b == *right_b && *left_c > *right_c) {
					//printf("      SWAP\n");
					int temp = *left_b;
					*left_b = *right_b;
					*right_b = temp;

					temp = *left_c;
					*left_c = *right_c;
					*right_c = temp;
					/*
				} else if (*left_b == *right_b && *left_c > *right_c) {
					int temp = *left_c;
					*left_c = *right_c;
					*right_c = temp;
					*/
				} else {
					//printf("\n");
					break;
				}
				//printf("After     : %i %i ... %i %i\n", *left_b, *right_b, *left_c, *right_c);
			}
		}
	}

	/*
	printf("INSIDE:\n");
	for (int i = 0; i < n * n; i++) {
		if (buckets_a[i] != 0)
			printf("%i ", buckets_a[i]);
	}
	printf("\n\n");

	for (int i = 0; i < n * n; i++) {
		if (buckets_b[i] != 0)
			printf("%i ", buckets_b[i]);
	}
	printf("\n\n");

	for (int i = 0; i < n * n; i++) {
		if (buckets_c[i] != 0)
			printf("%i ", buckets_c[i]);
	}
	printf("\n\n");
	*/


	int j = n - 1;
	j = 0;
	for (int i = 0; i < n * n; i++) {
		if (buckets_c[i] != 0) {
			dest[j++] = buckets_c[i];
			//printf("\t\tORDER: %i\n", buckets_c[i]);
		}
	}



	free(bucket_i);
	free(buckets_a);
	free(buckets_b);
	free(buckets_c);
}

/* Calculates the meet of two canonical factors A and B given as descending
 * cycle decomposition tables. */
void meet(uint8_t *dest, uint8_t *desc_a, uint8_t *desc_b, int n) {
	uint8_t *m = malloc(n * sizeof(uint8_t));
	for (int i = 0; i < n; i++)
		m[i] = n - i;

	sort_triples(m, desc_a, desc_b, m, n);

	/*
	printf("MEET:\n");
	for (int i = 0; i < n; i++)
		printf("%i ", m[i]);
	printf("\n");
	*/

	int j = m[n - 1];
	dest[j - 1] = j;

	/*
	for (int i = 0; i < n; i++)
		printf("%i ", m[i]);
	printf("\n");
	*/

	for (int i = n - 2; i >= 0; i--) {
		if (desc_a[j - 1] != desc_a[m[i] - 1]
				|| desc_b[j - 1] != desc_b[m[i] - 1]) {
			j = m[i];
		}
		dest[m[i] - 1] = j;
	}

	/*
	for (int i = 0; i < n; i++)
		printf("%i ", dest[i]);
	printf("\n");
	*/

	free(m);
}

void left_canonical_form(int8_t **dest, int *dest_len, int8_t *braid, int len, int n)
{
#if 0
	
	int as[8] = {2, 6, 3, 4, 1, 5, 7, 4};
	int bs[8] = {4, 2, 4, 4, 5, 7, 1, 2};
	int cs[8] = {8, 7, 6, 5, 4, 3, 2, 1};

	/*
	 * 712, 627, 573, 445, 421, 346, 248, 154
	 *
	 * 2, 7, 3, 5, 1, 6, 8, 4
	 */

	int result[8] = {0};
	/*
	sort_triples(result, as, bs, cs, 8);

	for (int i = 0; i < 8; i++) {
		printf("%i ", result[i]);
	}
	printf("\n");
	printf("\n");
	printf("\n");
	*/




	int a[8] = {5, 3, 3, 5, 5, 6, 7, 8};
	int b[8] = {4, 4, 3, 4, 5, 6, 7, 8};
	// cs[8] = {8, 7, 6, 5, 4, 3, 2, 1};
	meet(result, a, b, 8);


	printf("\n");
	printf("\n");
	for (int i = 0; i < 8; i++) {
		printf("%i ", result[i]);
	}
	printf("\n");




	printf("test\n");
	int pe[8] = {5, 3, 2, 1, 4, 6, 7, 8};
	int de[8] = {0};
	perm_to_desc(de, pe, 8);
	for (int i = 0; i < 8; i++) {
		printf("%i ", de[i]);
	}
	printf("\n");

	return;
#endif




	uint8_t trivial[8]       = {1, 2, 3, 4, 5, 6, 7, 8}; // Trivial permutation.
	//int delta[8]         = {8, 7, 6, 5, 4, 3, 2, 1}; // Fundamental braid...
	//int delta_inverse[8] = {8, 7, 6, 5, 4, 3, 2, 1}; // ... and its inverse.
	uint8_t delta[8]         = {8, 1, 2, 3, 4, 5, 6, 7}; // Fundamental braid...
	uint8_t delta_inverse[8] = {2, 3, 4, 5, 6, 7, 8, 1}; // ... and its inverse.

	/*
	int abc[8] = {0};
	perm_to_desc(abc, delta_inverse, 8);
	for (int i = 0; i < n; i++)
		printf("%i ", abc[i]);
	printf("\n");
	*/





	/*
	int perm[8] = {2, 4, 1, 3, 5, 6, 7, 8}; // 3 1 2
	inverse(perm, perm, 8);
	inverse(perm, perm, 8);

	for (int i = 0; i < 8; i++)
		printf("%i ", perm[i]);
	printf("\n");


	// tau
	multiply(perm, delta_inverse, perm, 8);
	multiply(perm, perm, delta, 8);
	// tau inverse
	multiply(perm, delta, perm, 8);
	multiply(perm, perm, delta_inverse, 8);

	for (int i = 0; i < 8; i++)
		printf("%i ", perm[i]);
	printf("\n");

	// negative substitution
	int neg[8] = {1, 2, 4, 3, 5, 6, 7, 8};
	multiply(neg, delta, neg, 8);
	multiply(neg, delta_inverse, neg, 8);

	for (int i = 0; i < 8; i++)
		printf("%i ", neg[i]);
	printf("\n");

	printf("\n");
	*/





	/*
	int testbraid[3] = {3, 4, -3};
	braid = testbraid;
	length = 3;
	*/

	uint8_t *band_braid = malloc(len * n * sizeof(uint8_t));
	int delta_power = 0;

	int *factors_delta_powers = malloc(len * sizeof(int));

	for (int i = 0; i < len; i++) {
		int abs_braid = abs(braid[i]);
		uint8_t *factor = band_braid + (i * n);
		// Convert into permutation table.
		for (int j = 0; j < n; j++) {
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

		for (int j = 0; j < n; j++)
			printf("%i ", factor[j]);
		printf(" (dp: %i)\n", factors_delta_powers[i]);
	}

	printf("\n");

	for (int i = len - 1; i > 0; i--) {
		uint8_t *right = band_braid + (i * n);
		uint8_t *left = right - n;

		int right_delta_power = factors_delta_powers[i];
		if (right_delta_power < 0) {
			// Apply index-shifting automorphism on left factor.
			for (int j = 0; j < abs(right_delta_power); j++) {
				multiply(left, delta, left, n);
				multiply(left, left, delta_inverse, n);
			}
			factors_delta_powers[i - 1] += right_delta_power;
		}

		for (int j = 0; j < n; j++)
			printf("%i ", right[j]);
		printf(" (dp: %i)\n", factors_delta_powers[i]);
	}
	delta_power += factors_delta_powers[0];
	for (int j = 0; j < n; j++)
		printf("%i ", band_braid[j]);
	printf(" (dp: %i)\n", factors_delta_powers[0]);

	printf("fundamental power: %i\n\n", delta_power);

	free(factors_delta_powers);

	/* next: Algorithm for converting permutation tables to the other thing.
	 * implement meet algorithm for the band generators.
	 */




	// Convert permutation tables into descending cycle decomposition tables.
	/*
	length = 1;
	int b[8] = {1, 2, 4, 3, 5, 6, 7, 8};
	band_braid = delta;
	multiply(band_braid, band_braid, b, 8);



	int *band_braid_desc = malloc(length * n * sizeof(int));
	for (int i = 0; i < length; i++)
		perm_to_desc(band_braid_desc + (i * n), band_braid + (i * n), n);

	for (int i = 0; i < n; i++)
		printf("%i ", band_braid_desc[i]);
	printf("\n");

	free(band_braid_desc);
	*/




	/*
	// test
	printf("\n\n\nFundamental Power: %i\n", delta_power);
	printf("Length: %i\n", length);
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < n; j++) {
			printf("%i ", band_braid[i * n + j]);
		}
		printf("\n");
	}
	printf("\n");
	*/


	

	// Convert braid into left canonical form.
	int i = 0;
	while (i < len) {
		int t = len;
		for (int j = len - 1; j > i; j--) {
			uint8_t *right = band_braid + (j * n);
			uint8_t *left = right - n;

			uint8_t *left_inverse = malloc(n * sizeof(uint8_t));
			inverse(left_inverse, left, n);

			uint8_t *left_star = malloc(n * sizeof(uint8_t));
			multiply(left_star, left_inverse, delta, n);

			// Convert operands to descending cycle decomposition tables.
			uint8_t *left_star_desc = malloc(n * sizeof(uint8_t));
			perm_to_desc(left_star_desc, left_star, n);
			uint8_t *right_desc = malloc(n * sizeof(uint8_t));
			perm_to_desc(right_desc, right, n);

			uint8_t *lr_meet_desc = malloc(n * sizeof(uint8_t));
			meet(lr_meet_desc, left_star_desc, right_desc, n);

			// Convert meet to permutation table.
			uint8_t *lr_meet = malloc(n * sizeof(uint8_t));
			desc_to_perm(lr_meet, lr_meet_desc, n);

			/*
			printf("ALGORITHM:\n");
			for (int i = 0; i < n; i++)
				printf("%i ", lr_meet[i]);
			printf("\n");
			for (int i = 0; i < n; i++)
				printf("%i ", left_star[i]);
			printf("\n");
			for (int i = 0; i < n; i++)
				printf("%i ", right_desc[i]);
			printf("\n");
			*/






			// --->
			/*
			bool nontrivial = false;
			for (int i = 0; i < n; i++) {
				if (lr_meet[i] != i) {
					nontrivial = true;
					break;
				}
			}
			*/
			if (!equal(lr_meet, trivial, n)) {
				t = j - 1;
				multiply(left, left, lr_meet, n);
				inverse(lr_meet, lr_meet, n);
				multiply(right, lr_meet, right, n);
			}



			free(left_star_desc);
			free(right_desc);
			free(lr_meet_desc);
			free(lr_meet);
			free(left_star);
			free(left_inverse);
		}
		i = t + 1;
	}



	// test
	printf("\nFundamental Power Before: %i\n", delta_power);
	printf("Length Before: %i\n", len);
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < n; j++) {
			printf("%i ", band_braid[i * n + j]);
		}
		printf("\n");
	}
	printf("\n");

	uint8_t *band_braid_reduced = band_braid;

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

	for (; len > 0
			&& equal(band_braid_reduced + n * (len - 1), trivial, n);
			len--) {
	}


	// test
	printf("\n\n\nFundamental Power After: %i\n", delta_power);
	printf("Length After: %i\n", len);
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < n; j++) {
			printf("%i ", band_braid_reduced[i * n + j]);
		}
		printf("\n");
	}
	printf("\n");

	// ADD FUNDAMENTAL POWERS TO BEGINNING!



	// Convert band braid into descending cycle decomposition tables.
	uint8_t *band_braid_desc = malloc(len * n * sizeof(uint8_t));
	for (int i = 0; i < len; i++)
		perm_to_desc(band_braid_desc + i * n, band_braid_reduced + i * n, n);

	printf("Descending cycles:\n");
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < n; j++) {
			printf("%i ", band_braid_desc[i * n + j]);
		}
		printf("\n");
	}
	printf("\n");


	free(band_braid);






	// Convert band braid into Artin generators.
	//
	// How big can it get??????
	int8_t *artins = malloc(len * n * n * sizeof(int8_t));
	int mark = 0;

	// Add fundamental braids to beginning.
	int8_t delta_artins[7] = {7, 6, 5, 4, 3, 2, 1};
	int8_t delta_inverse_artins[7] = {-1, -2, -3, -4, -5, -6, -7};
	for (int i = 0; i < abs(delta_power); i++) {
		if (delta_power > 0)
			memcpy(artins + mark, delta_artins, (n - 1) * sizeof(int8_t));
		else
			memcpy(artins + mark, delta_inverse_artins, (n - 1) * sizeof(int8_t));
		mark += n - 1;
	}
	
	for (int i = 0; i < len; i++) {
		uint8_t *factor = band_braid_desc + (i * n);
		bool *done = calloc(n, sizeof(bool));
		for (int j = n - 1; j >= 0; j--) {
			if (factor[j] !=  j + 1 || done[factor[j] - 1])
				continue;
			
			int t = j;
			int s = t - 1;
			for (; s >= 0; s--) {
				//printf("%i %i\n", factor[s], factor[t]);
				if (factor[s] != factor[t])
					continue;
				for (int k = 0; k < t - s; k++)
					artins[mark++] = t - k;
				for (int k = 0; k < t - s - 1; k++)
					artins[mark++] = -1 * (s + k + 2);
				t = s;
			}

			done[factor[t] - 1] = true;
		}
		free(done);
	}
	
	printf("\n\nArtin generators:\n");
	for (int i = 0; i < mark; i++) {
		printf("%i ", artins[i]);
	}
	printf("\n");



	free(band_braid_desc);

	//*artins_size = mark;
	//return artins;

	*dest = artins;
	*dest_len = mark;
}

#endif

