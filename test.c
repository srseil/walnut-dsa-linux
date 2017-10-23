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

//#include "galois.h"
#include "bkl.c"
#include "dehornoy.c"
#include "debug.c"

#include "dehornoy/braid.h"
#include "dehornoy/word.c"

#define PURE_BRAID_SIZE(x, n) (2 * ((n) - 1 - (x)) + 2)

/* Note that both the artin generators and the numbers in the permutations
 * start at 1 and not at 0!
 */

struct cb_pair {
	int n;				// Number of strands of the braid group dealt with.
	int q;				// Size of the galois field operated in.
	int *matrix;		// Pointer to an n*n matrix over the galois field.
	int *permutation;	// Pointer to an array of numbers modulo n.
};


//void emult(struct cb_pair *result, struct cb_pair m, int braid, int *t_values);
void emult(int *matrix, int *perm, int artin, int *t_values, int q, int n);
void gen_cloaking_elem(int **dest, int *dest_size,
	int *perm, int a, int b, int n, int bits);
void gen_pub_key(int **dest_matrix, int **dest_perm,
		int *priv_key, int priv_key_len, int *t_values, int q, int n);
int *gen_priv_key(int a, int b, int n, int bits);
void gen_priv_key2(int **dest, int *dest_size, int n, int bits);
void gen_pure_braid(int *dest, int x, int n);
void gen_pure_braids(int **dest, int *dest_size, int min_len, int n);
void encode(int **dest, int *dest_size, char *message, int size, int n);

void generate_signature(int **dest, int *dest_len, char *hash, int hash_len,
		int *priv_key, int priv_key_len, int *pub_key_matrix, int *pub_key_perm,
		int *t_values, int bits);
int verify_signature(int *sig, int sig_len, char *hash, int hash_len,
		int *pub_key_matrix, int *pub_key_perm, int *t_values, int q, int n);

void example();
void cloaking_example();


uint8_t galois_mult(uint8_t a, uint8_t b) {
	uint8_t r = 0, t;
	while (a != 0) {
		if ((a & 1) != 0)
			r ^= b;
		t = b & 0x10;
		b <<= 1;
		if (t != 0)
			b ^= 0x25;
		a >>= 1;
	}
	return r;
}


uint8_t galois_inverse(uint8_t x) {
	for (int i = 0; i < 32; i++) {
		if (galois_mult(x, i) == 1u)
			return i;
	}
	return 0;
}




void dehornoy_old(int **dest, int *dest_len, int *braid, int len) {
	// First, find fitting pairs of generators in braid to reduce.
	// If we delete something (esp. freely reduce!), replace with zeroes

	// Make sure that this braid has enough space for the reduction expansions!
	/*
	int braid_inline[] = {
		//1, 0, 0, 2, 0, 0, 0, -2, 3, 3, 0, -3, 4, 5, -5
		//
		-1, -2, 1, 3, -2, -3, -2, 1, -3, 2, 1, 1
		//0, -1, 0, -2, 0, 0, 1, 3, -2, 0, -3, -2, 1, -3, 2, 1, 1

		//2, -1, -2, 3, -2, -3, -2, 1, -3, 2, 1, 1
		//2, 0, -1, -2, 0, 0, 3, -2, 0, -3, -2, 1, 0, -3, 2, 1, 1
		
		//2, 2, -1, 3, 2, -1, -2, -3, -1, 3, -2, -3, 1, 1
		//2, 2, 0, 0, -1, 3, 0, 2, -1, -2, 0, -3, -1, 0, 3, 0, -2, -3, 1, 1

		//2, 2, -1, 3, 2, -1, -2, -3, -1, -2, -3, 2, 1, 1

		//2, 2, -1, 3, 2, -1, 2, -1, -2, -3 // no handle
	};
	len = sizeof(braid_inline) / sizeof(int);
	braid = calloc(2 * len, sizeof(int));
	memcpy(braid, braid_inline, sizeof(braid_inline));
	//braid = braid_inline;
	*/


	int orig = len;
	int *temp = braid;
	braid = calloc(10 * len, sizeof(int));
	memcpy(braid, temp, len * sizeof(int));


	//printf("Braid: ");
	//print_braid(braid, len);

	// Maybe do free reduction once before this?
	for (int i = 0; i < len - 1; i++) {
		if (braid[i] == -1 * braid[i + 1]) {
			braid[i] = 0;
			braid[i + 1] = 0;
			i++;
		}
	}

	int cnt = 0;

	// Handle reduction
	while (true) {

	// We have to reduce the nested handle first, if there is any...
	// Find handle.
	
	int search_start = 0, search_end = len;
	int start = -1, end = -1;
	int handle_found = 1;
	while (handle_found) {
		handle_found = 0;
		for (int i = search_start; i < search_end - 2; i++) {
			if (braid[i] == 0)
				continue;
			for (int j = i + 1; j < search_end; j++) {
				if (braid[j] == 0) {
					continue;
				} else if (braid[j] == braid[i]
						|| abs(braid[j]) == abs(braid[i]) - 1) {
					break;
				} else if (braid[j] == -1 * braid[i] && (j - i) > 1) {
					start = i;
					end = j;
					search_start = i + 1;
					search_end = j + 1;
					handle_found = 1;
					goto end;
				}
			}
		} end:;
	}

	if (start == -1 || end == -1) {
		printf("No handle found.\n");
		break;
	}

	/*
	printf("Handle: (%i) - (%i) (pos %i-%i)\n", braid[start], braid[end],
		start, end);
		*/

	// 2 * too short for handle?
	int *handle = calloc(10 * (end - start + 1), sizeof(int));
	if (handle == NULL)
		puts("handle is NULL");

	int start_sig = sig(braid[start]);
	int main_gen = abs(braid[start]);
	//printf("start_sig: %i, main_gen: %i\n", start_sig, main_gen);
	int off = 0;

	// Reduce handle.
	int j = 1;
	for (int i = 1; i < end - start; i++) {
		if (braid[start + i] == 0)
			continue;
		if (abs(braid[start + i]) == main_gen + 1) {
			handle[j + 2 * off - 1] = -1 * start_sig * (main_gen + 1);
			handle[j + 2 * off] = sig(braid[start + i]) * main_gen;
			handle[j + 2 * off + 1] = start_sig * (main_gen + 1);
			off++;
		} else {
			handle[j + 2 * off - 1] = braid[start + i];
		}
		j++;
	}

	//printf("Handle: ");
	//print_braid(handle, 2 * (end - start + 1)); // FALSCHE LÄNGE


	//int handle_len = end - start + 2 * off - 1;
	int handle_len = j + 2 * off - 1;
	//printf("handle_len: %i\n", handle_len);
	printf("handle_len: %i\n", j + 2 * off - 1);

	for (int i = start; i <= end; i++)
		braid[i] = 0;

	if (start + handle_len > end) {
		memmove(braid + start + handle_len, braid + end + 1,
				(len - end - 1) * sizeof(int));
		printf("len before: %i, ", len);
		len += start + handle_len - end - 1;
		printf("after: %i      (orig: %i)\n", len, orig);
	}

	// memcpy?
	for (int i = 0; i < handle_len; i++)
		braid[start + i] = handle[i];

	free(handle);

	//print_braid(braid, len);


	// Free-reduce braid.
	for (int i = 0; i < len - 1; i++) {
		if (braid[i] == 0)
			continue;
		int j = i + 1;
		for (; j < len && braid[j] == 0; j++);

		if (j < len && braid[i] == -1 * braid[j]) {
			braid[i] = 0;
			braid[j] = 0;
			i = j;
		}
	}

	/*
	if (++cnt == 20) {
		printf("\nBEFORE:\n");
		print_braid(braid, len);
	}

	5,6
	2,4
	0,2
	*/

	// Remove zeroes to shrink braid.
	int i = 0;
	for (; i < len && braid[i] != 0; i++);

	for (int j = i + 1; j < len; j++) {
		if (braid[j] != 0)
			braid[i++] = braid[j];
	}
	
	len = i;

	/*
	if (cnt == 20) {
		printf("\nAFTER:\n");
		print_braid(braid, len);
		cnt = 0;
	}
	*/

	if (++cnt == 20) {
		print_braid(braid, len);
		cnt = 0;
	}
	
	}

	int *result = malloc(len * sizeof(int));
	int mark = 0;
	for (int i = 0; i < len; i++) {
		if (braid[i] != 0)
			result[mark++] = braid[i];
	}

	*dest = result;
	*dest_len = mark;
}


int main(int argc, char *argv[]) {

#if 0
	//CREATE VALUE TABLE FOR GALOIS INVERSE IN G_256
	for (int i = 0; i < 256; i++) 
		printf("%i, ", galois_inverse(i));
	puts("\n");
	return 0;
#endif


#if 0
	/*int schmeckle[] = {
		-1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, -1, -2, -3, -4, -5, -6, -7, 7, 6, 5, 4, 3, 1, 7, 6, 5, 4, 3, 2, 1, -2, -3, 2, 7, 6, 5, 4, 3, 2, 1, -2, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, 4, 3, 2, -3, 7, 6, 5, 4, -5, -6, 3, 2, -3, 1, 7, 6, 5, 4, 3, 2, -3, -4, -5, -6, -7, 1, 6, 5, 4, 3, -4, -5, -6, 5, 4, -5, 7, 6, 5, -6, -7, 6, 3, 2, 1, 7, 6, 5, 4, -5, -6, 3, 2, 1, 5, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, -7, 6, 4, 3, 2, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, -7, 5, 4, 3, -4, -5, 6, 3, 7, 5, 4, -5, 7, 7, 6, 5, 4, -5, -6, -7, 7, 6, 5, 4, 3, 2, -3, -4, -5, -6, -7, 6, 5, 4, -5, -6, 7, 5, 4, -5, 3, 2, 1, -2, -3, 2, 7, 6, -7, 5, 4, 3, 2, 1, -2, -3, -4, -5, 4, 3, 6, 5, -6, 4, 3, 2, -3, 1, 7, 6, 5, 4, -5, -6, 3, 2, -3, 5, 6, 5, 4, 3, -4, -5, 2, 1, -2, 7, 6, 5, -6, 3, 2, -3, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, 5, 4, 3, 2, -3, -4, -5, 7, 6, 5, 4, 3, 2, -3, -4, -5, -6, -7, 1, 6, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, 2, 7, 6, 5, -6, 3, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, 5, 3, 2, -3, 7, 6, 4, 3, 2, 1, -2, -3, 2, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, 5, 4, 3, -4, -5, 2, 4, 7, 6, 5, 4, 3, 2, -3, -4, -5, -6, -7, 1, 6, 5, 4, 3, -4, -5, 7, 6, 5, 4, -5, -6, -7, 3, 2, -3, 1, 6, 7, 4, 3, 2, 1, -2, -3, -4, 3, 2, -3, 7, 6, 5, 4, 3, 2, -3, -4, -5, -6, -7, 1, 4, 4, 3, -4, 1, 7, 6, 5, 4, -5, -6, -7, 5, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, -7, 6, 5, 4, 3, 2, -3, -4, -5, -6, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, -7, 3, 2, -3, 7, 6, 5, 4, -5, -6, -7, 3, 2, 1, -2, -3, 2, 7, 6, 5, 4, -5, -6, -7, 3, 2, 1, -2, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, -7, 6, 5, 4, 3, -4, -5, -6, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, -7, 6, 5, 4, 3, -4, -5, -6, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, -7, 6, 5, 4, 3, 2, -3, -4, -5, -6, 4, 3, -4, 7, 6, 5, 4, 3, -4, -5, -6, -7, 5, 5, 5, 7, 6, -7, 7, 6, -7, 7, 6, -7, 7, 6, -7, 7, 6, -7, 7, 6, -7, 7, 6, -7, 7, 6, 5, -6, -7, 6, 7, 6, 5, -6, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, -7, 6, 5, 4, -5, -6, 5, 6, 2, 1, -2, 7, 6, 5, 4, 3, 2, -3, -4, -5, -6, -7, 6, 5, 4, -5, -6, 6, 5, 4, -5, -6, 2, 1, -2, 7, 6, 5, 4, 3, 2, -3, -4, -5, -6, -7, 6, 5, 4, -5, -6, 5, 4, -5, 3, 2, 1, -2, -3, 2, 6, 4, 3, -4, 2, 1, 6, 5, 4, 3, -4, 1, 7, 6, 5, 4, 3, 2, -3, -4, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, 4, 3, 2, -3, -4, 7, 5, 4, 3, 2, 1, -2, -3, -4, 2, 7, 5, 2, 1, 7, 4, 3, -4, 1, 7, 5, 2, 7, 5, 4, 3, -4, -5, 2, 7, 6, 5, 4, 3, 2, -3, -4, -5, -6, -7, 6, 4, 3, -4, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, 5, 4, 3, 2, -3, -4, 3, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, -7, 6, 5, 4, 2, 6, 5, 4, 3, -4, 2, 1, 7, 6, 5, 4, -5, 3, 7, 6, -7, 5, 4, 3, -4, -5, 5, 4, 3, -4, -5, 6, 3, 6, 5, 4, -5, 4, 7, 6, 5, -6, -7, 7, 6, 5, -6, -7
	};
	*/

	int schmeckle[] = {
		/*
		7, -6, 5, -6, 5, -6, 5, -4, 5, -7, 6, 7, 5, -4, -4, -5, 6, -4, -3, -7, 6, 4, 5, 6, 5, -3, -3, 2, 1, -6, 7, -6, -7, -6, 5, 6, -4, 3, -2, 1, 4, -3, 2, 5, -4, -5, 3, -4, 5, 2, -3, 4, 1, -2, 1, 6, 5, -6, 4, 3, 2, 2, 1, 4, 3, 2, 1, 2, 1, 5, 4, 3, 2, 2, 3, 2, 1, 4, 3, 2, 7, 6, 5, -6, -7, 4, -5, -6, 3, -4, -5, 2, -3, -4, 1, -2, -3, -2, -2, -2, -2, -2, -2, -2, -2, -3, -4, -5, -6, -7, 1, 1, -3, -4, -5, -6, -4, -4, -5, -6, -5, -5, -3, -3, 1, 1, -3, -4, -5, -4, -4, -5, -4, -4, 5, -6, -4, 2, 1, 3, 2, 1, 1, 4, -5, 3, -4, 3, 2, 2, 3, 3, 2, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -3, -2, -2, -3, -3, -2, -2, 1, -3, 2, 6, -5, 7, -6, -7, -5, -4, 3, 4, 5, 2, 3, 4, 1, 2, 3, 2, 1, 2, -5, -4, 5, -6, 3, -4, -5, 3, -4, 5, 6, 2, 1, 1, 1, -2, -3, 1, 1, -2, -3, -4, -5, -6, 1, 1, 1, 1, 1, 1, 1, -2, -3, -4, -5, 6, 6, 6, -3, -4, 5, 5, 6, 7, 5, 6, 6, 7, 5, 6, 6, 6, 5, 5, -2, -3, -2, -2, -3, -4, -5, -6, -7, 1, -2, -3, -4, -5, -6, -2, -3, -4, -5, -4, -3, -3, -4, -2, -3, 4, -3, -2, -6, -6, -6, 7, 7, 7, 7, 7, 7, 7, -6, -7, 1, 1, 2, 3, -4, 2, -3, 1, -2, 1, 5, 4, 3, 2, -3, -4, -5, -3, -4, 5, 1, -2, 3, 6, 5, 4, -5, -6, 3, -4, -4, 7, -6, 1, 2, 5, 3, 2, -4, 3, 4, 1, 6, -7, 5, -6, 4, 3, -2, 1, 4, -5, -5, -6, -6, 7, -6, -7, 3, 3, -4, 5, 3, 5, -6, 4, -5, 4, 7, 6, 5
			*/

			7, -6, 5, -6, 5, -6, 5, -4, 5, -7, 6, 7, 5, -4, -4, -5, 6, -4, -3, -7, 6, 4, 5, 6, 5, -3, -3, 2, 1, -6, 7, -6, -7, -6, 5, 6, -4, 3, -2, 1, 4, -3, 2, 5, -4, -5, 3, -4, 5, 2, -3, 4, 1, -2, 1, 6, 5, -6, 4, 3, 2, 2, 1, 4, 3, 2, 1, 2, 1, 5, 4, 3, 2, 2, 3, 2, 1, 4, 3, 2, 7, 6, 5, -6, -7, 4, -5, -6, 3, -4, -5, 2, -3, -4, 1, -2, -3, -2, -2, -2, -2, -2, -2, -2, -2, -3, -4, -5, -6, -7, 1, 1, -3, -4, -5, -6, -4, -4, -5, -6, -5, -5, -3, -3, 1, 1, -3, -4, -5, -4, -4, -5, -4, -4, 5, -6, -4, 2, 1, 3, 2, 1, 1, 4, -5, 3, -4, 3, 2, 2, 3, 3, 2, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -3, -2, -2, -3, -3, -2, -2, 1, -3, 2, 6, -5, 7, -6, -7, -5, -4, 3, 4, 5, 2, 3, 4, 1, 2, 3, 2, 1, 2, -5, -4, 5, -6, 3, -4, -5, 3, -4, 5, 6, 2, 1, 1, 1, -2, -3, 1, 1, -2, -3, -4, -5, -6, 1, 1, 1, 1, 1, 1, 1, -2, -3, -4, -5, 6, 6, 6, -3, -4, 5, 5, 6, 7, 5, 6, 6, 7, 5, 6, 6, 6, 5, 5, -2, -3, -2, -2, -3, -4, -5, -6, -7, 1, -2, -3, -4, -5, -6, -2, -3, -4, -5, -4, -3, -3, -4, -2, -3, 4, -3, -2, -6, -6, -6, 7, 7, 7, 7, 7, 7, 7, -6, -7, 1, 1, 2, 3, -4, 2, -3, 1, -2, 1, -2, -7, 5, 4, 3, -4, -5, 2, -3, -4, -4, -5, 1, -2, -3, -4, -4, -5, -5, -6, -3, -2, 4, 7, -5, 2, 1, -2, -2, -3, -4, -5, -6, 5, 4, 3, 2, -3, -4, 3, 7, 6, 5, 4, 3, 2, 1, -2, -3, -4, -5, -6, -7, 6, 5, 4, 2, 6, 5, 4, 3, -4, 2, 1, 7, 6, 5, 4, -5, 3, 7, 6, -7, 5, 4, 3, 3, -4, -5, 6, 3, 6, 5, 4, -5, 4, 7, 6, 5
	};

	int schlen = sizeof(schmeckle) / sizeof(int);
	printf("len: %i\n", schlen);

	int per[] = {1, 2, 3, 4, 5, 6, 7, 8};
	for (int i = 0; i < schlen; i++) {
		int k = abs(schmeckle[i]);
		int temp = per[k - 1];
		per[k - 1] = per[k];
		per[k] = temp;
	}

	for (int i = 0; i < 8; i++)
		printf("%i ", per[i]);
	printf("\n");

	return 0;
#endif





	//test_shit();
	//dehornoy(NULL, NULL, NULL, 0);
	//return 0;

#if 0
	int braid_inline[] = {
		//-1, -2, 1, 3, -2, -3, -2, 1, -3, 2, 1, 1
		4, 5, -2, 7, -5, 7, -4, 7, -6, -6, 7, 1, 1, 2, 4, -5, 7, 7, 6, 6, 1, -5,
		7, 5, 3, 1, 5, 2, 5, 7, 5, -4, 5, -2, -4, 2, 4, -1, -1, -6, -4, -1, -2,
		7, -5, 1, -2, 4, 6, -1, -5, 6, 2, -1, 6, -7, -2, -3, -2, -5, -4, -7, -1,
		2, 2, 4, -6, 2, -4, 1, 1, 7, 6, -5, 3, 7, 2, -6, -5, 1, 5, 4, 3, 6, -7,
		-2, 4, 4, -6, 3, -5, 6, 3, 7, -5, -4, -7, -1, -6, -6, 7, -6, -5, -7, 5,
		3, 5, 5, 1, 4, -6, -5, 4, 1
	};
	int *braid_res, braid_res_len;
	dehornoy(&braid_res, &braid_res_len, braid_inline, sizeof(braid_inline) / sizeof(int));
	print_braid(braid_res, braid_res_len);

	int per[] = {1, 2, 3, 4, 5, 6, 7, 8};
	//for (int i = 0; i < sizeof(braid_inline) / sizeof(int); i++) {
		//int k = abs(braid_inline[i]);
	for (int i = 0; i < braid_res_len; i++) {
		int k = abs(braid_res[i]);
		int temp = per[k - 1];
		per[k - 1] = per[k];
		per[k] = temp;
	}

	for (int i = 0; i < 8; i++)
		printf("%i ", per[i]);
	printf("\n");

	// 36152784
	return 0;
#endif	

#if 0
	// Galois test.
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

	//int tt_values[8] = {1, 1, 20, 14, 5, 9, 10, 7};
	int tt_values[8] = {14, 1, 20, 1, 5, 9, 10, 7};
	int tbraid[] = {
		4, 5, -2, 7, -5, 7, -4, 7, -6, -6, 7, 1, 1, 2, 4, -5, 7, 7, 6, 6, 1, -5,
		7, 5, 3, 1, 5, 2, 5, 7, 5, -4, 5, -2, -4, 2, 4, -1, -1, -6, -4, -1, -2,
		7, -5, 1, -2, 4, 6, -1, -5, 6, 2, -1, 6, -7, -2, -3, -2, -5, -4, -7, -1,
		2, 2, 4, -6, 2, -4, 1, 1, 7, 6, -5, 3, 7, 2, -6, -5, 1, 5, 4, 3, 6, -7,
		-2, 4, 4, -6, 3, -5, 6, 3, 7, -5, -4, -7, -1, -6, -6, 7, -6, -5, -7, 5,
		3, 5, 5, 1, 4, -6, -5, 4, 1
	};
	int tbraid_size = sizeof(tbraid) / sizeof(int);
	/*
	int *tbraid;
	int tbraid_size;
	gen_priv_key2(&tbraid, &tbraid_size, 8, 64);
	*/


	int *pubk_matrix, *pubk_perm;
	gen_pub_key(&pubk_matrix, &pubk_perm, tbraid, tbraid_size,
	//gen_pub_key(&pubk_matrix, &pubk_perm, tbraid, sizeof(tbraid) / sizeof(int),
		tt_values, 32, 8);

	int *signature, signature_len;
	char hash[4] = {0x11, 0x22, 0x33, 0x44};
	// Hash "length"? More like size, since it means the bytes?
	generate_signature(&signature, &signature_len, hash, 4, tbraid,
		//sizeof(tbraid) / sizeof(int), pubk_matrix, pubk_perm, tt_values, 64);
		tbraid_size, pubk_matrix, pubk_perm, tt_values, 64);

	verify_signature(signature, signature_len, hash, 4, pubk_matrix, pubk_perm, tt_values, 32, 8);

	return 0;






	int n = 8;
	int *identity = malloc(n * n * sizeof(int));
	for (int i = 0; i < n * n; i++) {
		identity[i] = (i % (n + 1) == 0) ? 1 : 0;
	}
	/*
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			printf("%i ", identity[i * n + j]);
		printf("\n");
	}
	printf("\n");
	*/

	int *id_perm = malloc(n * sizeof(int));
	for (int i = 0; i < n; i++)
		id_perm[i] = i + 1;


	for (int i = 0; i < sizeof(tbraid) / sizeof(int); i++) {
		emult(identity, id_perm, tbraid[i], tt_values, 32, n);
	}

	return 0;






	
	int *cloak = NULL;
	int cloak_size = 0;
	int cloak_perm[8] = {3, 2, 8, 5, 6, 7, 1, 4};
	gen_cloaking_elem(&cloak, &cloak_size, cloak_perm, 3, 6, 8, 32);
	print_braid(cloak, cloak_size);
	return 0;



	int *pri = NULL;
	int pri_len = 0;
	gen_priv_key2(&pri, &pri_len, 8, 64);
	print_braid(pri, pri_len);
	return 0;



	int *pures;
	int pures_len = 0;
	int min_len = 10;
	/*
	gen_pure_braids(&pures, &pures_len, min_len, 8);
	printf("min_len: %i, actual_len: %i\n", min_len, pures_len);
	print_braid(pures, pures_len);
	*/

	encode(&pures, &pures_len, "ab1d3fgAui89", 12 * sizeof(char), 8);
	print_braid(pures, pures_len);

	int purperm[] = {1, 2, 3, 4, 5, 6, 7, 8};
	for (int i = 0; i < pures_len; i++) {
		int k = abs(pures[i]);
		int temp = purperm[k - 1];
		purperm[k - 1] = purperm[k];
		purperm[k] = temp;
	}

	for (int i = 0; i < 8; i++)
		printf("%i ", purperm[i]);
	printf("\n");

	return 0;



	return 0;





	int key[] = {
		4, 5, -2, 7, -5, 7, -4, 7, -6, -6, 7, 1, 1, 2, 4, -5, 7, 7, 6, 6, 1, -5,
		7, 5, 3, 1, 5, 2, 5, 7, 5, -4, 5, -2, -4, 2, 4, -1, -1, -6, -4, -1, -2,
		7, -5, 1, -2, 4, 6, -1, -5, 6, 2, -1, 6, -7, -2, -3, -2, -5, -4, -7, -1,
		2, 2, 4, -6, 2, -4, 1, 1, 7, 6, -5, 3, 7, 2, -6, -5, 1, 5, 4, 3, 6, -7,
		-2, 4, 4, -6, 3, -5, 6, 3, 7, -5, -4, -7, -1, -6, -6, 7, -6, -5, -7, 5,
		3, 5, 5, 1, 4, -6, -5, 4, 1

		/*
		-4, 3, -6, 1, -4, -3, 1, -4, -5, 2, -1, -2, -2, -3, -1, 5, 5, -1, 4, 4,
		5, 4, 3, -4, -7, -3, 6, 5, -7, 4, 2, 1, -7, -5, 1, 2, 6, -4, -1, -2, -5,
		-1, 4, -3, 6, -3, 5, 1, -5, -2, 4, -6, -7, -1, -1, 3, -7, -4, -3, 2, -5,
		2, 5, -1, 6, 4, 2, 3, 4, -3, 2, 3, -4, 5, 6, 5, -4, 2, -1, 6, -7, -6
			*/

		/*
		2, 3, 4, 5, 6, -5, -4, 3, 2, 6, -4, 5, -6, -5, 4, 3, -4, -5, 6, 5, -4,
		-3, 1, -2, -3, -4, -5, 6, 7, -6, -5, -4, 3, -2, 1, 5, -6, -7, 6, 4, -3,
		-3, -4, -5, 2, -1, -1, -2, -6, -6 /*, 7, 7, 6, 6, 2, 1, 1, -2, 5, 4, 3, 3, 		-4, -6, 7, 6, -5, -1, 2, -3, 4, 5, 6, -7, -6, 5, 4, 3, 2, -1, 3, 4, -5,
		-6, 5, 4, -3, -4, 5, 6, -5, 4, -6, -2, -3, 4, 5, -6, -5, -4, -3, -2*/
	};

	int perm[] = {1, 2, 3, 4, 5, 6, 7, 8};
	for (int i = 0; i < sizeof(key) / sizeof(int); i++) {
		int k = abs(key[i]);
		int temp = perm[k - 1];
		perm[k - 1] = perm[k];
		perm[k] = temp;
	}

	for (int i = 0; i < 8; i++)
		printf("%i ", perm[i]);
	printf("\n");

	return 0;



	cloaking_example();
	return 0;

	//encode(int **dest, int *dest_size, char *message, int size, int n) {
	char message[] = {
		/*
		0x79, 0xb7, 0xac, 0x30, 0x39, 0x4a, 0xff, 0x82, 0x92, 0xed, 0x52, 0xf3,
		0x8d, 0x52, 0x0c, 0xf8 ,0x2b, 0x01, 0xc8, 0x00, 0x63, 0x07, 0x46, 0xef,
		0x61, 0x4b, 0x26, 0xf0, 0x45, 0xc1, 0x09, 0x5d
		*/
		0x2d, 0x18, 0x06, 0x35, 0x27, 0x39, 0x19, 0x91, 0x16, 0x7b, 0x12, 0x81,
		0x85, 0xde, 0xbe, 0x56, 0x7b, 0x3d, 0xd0, 0x1b, 0x01, 0xe2, 0xaf, 0x03,
		0xed, 0xce, 0x2d, 0xef, 0x04, 0x8e, 0x06, 0xcc
	};
	int *encoded = NULL;
	int enc_size = 0;
	encode(&encoded, &enc_size, message, 32, 8);


	for (int i = 0; i < enc_size; i++)
		printf("%i ", encoded[i]);
	printf("\n");


	/*
void gen_cloaking_elem(int **dest, int *dest_size,
		int a, int b, int c, int n, int bits) {

	int *clk, clk_size;
	gen_cloaking_elem(&clk, &clk_size, 3, 6, 
	*/

	return 0;



	// BKL test
	//int braid[] = {3, -4, 2, -6};

	int braid[] = {
		6,7,-6,1,-2,4,-5,-6,-5,4,-3,-2,3,-4,-3,-2,-4,-6,1,-5,-2,5,-2,3,4,7,-3,1,1,7,6,-4,2,5,-1,-5,3,-6,3,-4,1,5,2,1,4,-6,-2,-1,5,7,-1,-2,-4,7,-5,-6,3,7,4,-3,-4,-5,-4,-4,1,-5,-5,1,3,2,2,1,-2,5,4,-1,3,4,-1,6,-3,4,2,-3,4,5,4,3,-1,-2,3,-4,-5,-6,7,6,-5,4,3,2,1,-7,-5,-6,7,6,5,-4,-5,6,7,-6,5,4,2,-1,-1,-2,-7,6,5,4,3,2,-1,-1,-2,-3,-4,-5,-6,-7,-4,-4,-5,-5,-7,-7,6,5,-4,-4,-5,-6,3,-2,-2,-3,5,-4,-4,-5,-2,-2,5,4,-3,-3,-4,-5,7,6,5,4,-3,-3,-4,-5,-6,-7,6,-5,-5,-6,5,4,-3,-3,-4,-5,7,-6,-6,-7,-7,-7,-4,-4,-7,-7,2,2,7,7,4,4,7,7,7,6,6,-7,5,4,3,3,-4,-5,6,5,5,-6,7,6,5,4,3,3,-4,-5,-6,-7,5,4,3,3,-4,-5,2,2,5,4,4,-5,3,2,2,-3,6,5,4,4,-5,-6,7,7,5,5,4,4,7,6,5,4,3,2,1,1,-2,-3,-4,-5,-6,7,2,1,1,-2,-4,-5,6,-7,-6,5,4,-5,-6,-7,6,5,7,-1,-2,-3,-4,5,-6,-7,6,5,4,-3,2,1,-3,-4,-5,-4,3,-2,7,6,5,5,5,4,3,3,3,3,3,3,3,3,3,3,3,2,1,1,1,1,1,1,1,1,-2,-3,-4,5,5,5,-6,7,7,6,5,4,3,3,3,3,-4,5,-6,7,7,7,7,7,7,6,5,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,-4,5,5,5,-6,7,7,7,7,7,7,7,7,7,7,6,5,4,3,3,-4,5,5,4,3,2,1,1,1,1,1,1,-2,3,3,2,1,1,1,1,1,1,-2,3,3,3,3,3,3,3,3,3,3,3,-4,5,5,5,5,5,5,5,-6,7,7,7,7,7,7,6,5,5,5,5,5,5,5,5,5,4,3,3,3,3,-4,5,5,5,-6,7,7,7,7,7,7,7,7,7,7,7,7,6,5,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,1,1,-2,3,-4,-5,-6,7,7,7,7,7,7,6,5,4,3,2,1,1,-2,3,-4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,-6,7,7,7,7,7,7,7,7,6,5,4,3,2,1,1,-2,-3,-4,-5,-6,7,7,6,5,5,5,5,5,5,5,5,5,4,3,3,3,3,3,3,3,3,3,2,1,1,1,1,1,1,1,1,-2,-3,-4,5,5,5,5,5,5,5,5,5,5,4,3,3,3,3,3,3,3,3,-4,5,5,5,5,5,5,5,-6,7,7,7,7,7,7,7,7,6,5,4,3,2,1,1,1,1,1,1,1,1,1,1,1,1,-2,-3,-4,5,5,5,5,5,5,5,5,4,3,2,1,1,-2,-3,-4,5,5,5,5,4,3,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-2,-3,-4,-5,-6,-7,-4,3,-6,1,-4,-3,1,-4,-5,2,-1,-2,-2,-3,-1,5,5,-1,4,4,5,4,3,-4,-7,-3,6,5,-7,4,2,1,-7,-5,1,2,6,-4,-1,-2,-5,-1,4,-3,6,-3,5,1,-5,-2,4,-6,-7,-1,-1,3,-7,-4,-3,2,-5,2,5,-1,6,4,2,3,4,-3,2,3,-4,5,6,5,-4,2,-1,6,-7,-6,4,5,-6,-7,6,-5,4,5,-6,-7,6,5,1,2,-3,-4,-5,-6,-5,-4,-3,2,2,-3,4,5,-4,3,-2,-2,7,6,-5,-5,-6,-7,2,-1,-1,-2,6,5,4,-3,-3,-4,-5,-6,7,6,5,4,3,-2,-2,-3,-4,-5,-6,-7,6,5,4,3,2,-1,-1,-2,-3,-4,-5,-6,-7,-7,-3,-3,-6,-6,3,2,-1,-1,-2,-3,5,4,3,-2,-2,-3,-4,-5,-4,-4,6,-5,-5,-6,7,6,-5,-5,-6,-7,-7,-7,-6,5,4,3,-2,-2,-3,-4,-5,-6,7,6,5,4,3,2,-1,-1,-2,-3,-4,-5,-6,-7,5,4,3,2,-1,-1,-2,-3,-4,-5,2,2,5,4,3,2,1,1,-2,-3,-4,-5,7,6,5,4,3,2,1,1,-2,-3,-4,-5,-6,-7,6,5,4,3,2,2,-3,-4,-5,6,7,7,7,6,5,5,-6,-7,6,5,5,-6,4,4,5,4,3,2,2,-3,-4,-5,3,2,1,1,-2,-3,6,6,3,3,7,7,6,5,4,3,2,1,1,-2,-3,-4,-5,-6,7,6,5,4,3,2,2,-3,-4,-5,-6,-7,6,5,4,3,3,-4,-5,-6,2,1,1,-2,7,6,5,5,-6,-7,2,2,-3,4,-5,-4,3,-2,-2,3,4,5,6,5,4,3,-2,-1,-5,-6,7,6,-5,-4,5,-6,7,6,-5,-4
	};

	printf("size: %i\n", sizeof(braid) / sizeof(int));
	int braid_size = sizeof(braid) / sizeof(int);

	//int braid[] = {-2, -1, 2, -6};
	int out_size = 0;
	int *lcf = left_canonical_form(&out_size, braid, braid_size, 8);
	puts("\n");
	printf("BKL size: %i\n", out_size);
	puts("\n");


	/*
	// Dehornoy form with the fucked up library.
	LENGTH = out_size;
	memcpy(Word, lcf, out_size * sizeof(int));

	int *final = Reduction();

	//final[1] = -7;
	for (int i = 0; i < LENGTH; i++) {
		printf("%i ", final[i]);
	}
	printf("\n");

	printf("\n\n\nsize: %i -> %i\n", out_size, LENGTH);
	*/



	struct cb_pair m;
	m.n = 8;
	m.q = 32;

	int m_matrix[8*8] = {
		15, 30,  7, 18, 13, 20, 15, 31,
		10, 19, 13, 19,  6, 17, 11, 21,
		10, 14, 16, 19,  6, 17, 11, 21,
		17, 31, 21, 28, 15, 23,  2, 16,
		 9, 16, 10, 13, 12,  7, 31, 20,
		 9, 16, 10, 13, 21, 30, 31, 20,
		 0,  0,  2,  4, 23,  0,  8, 24,
		 0,  0,  0,  0,  0,  0,  0,  1
	};
	m.matrix = malloc(8 * 8 * sizeof(int));
	memcpy(m.matrix, &m_matrix, 8 * 8 * sizeof(int));

	int m_permutation[8] = {2, 4, 3, 8, 1, 6, 7, 5};
	m.permutation = malloc(8 * sizeof(int));
	memcpy(m.permutation, &m_permutation, 8 * sizeof(int));

	int tvals[8] = {28, 24, 1, 9, 26, 1, 18, 18};

	for (int i = 0; i < out_size; i++) {
		struct cb_pair result;
		result.n = 8;
		result.q = 32;
		int r_matrix[8*8] = {0};
		result.matrix = r_matrix;
		int r_permutation[8] = {0};
		result.permutation = r_permutation;

		//emult(&result, m, lcf[i], tvals);

		memcpy(m.matrix, result.matrix, 8 * 8 * sizeof(int));
		memcpy(m.permutation, result.permutation, 8 * sizeof(int));
	}

	free(m.matrix);
	free(m.permutation);


	return 0;
	// end

	//cloaking_example();
	//int *priv = gen_priv_key(3, 6, 8, 128);
	//example();

	/*
	char msg[] = {121u, 160u}; // 0111 1001 1010 0000
	int *enc, enc_size;
	encode(&enc, &enc_size, msg, 2, 8);

	for (int i = 0; i < enc_size; i++)
		printf("%i, ", enc[i]);
	printf("\n");
	*/

	return 0;
}

void gen_pub_key(int **dest_matrix, int **dest_perm,
		int *priv_key, int priv_key_len, int *t_values, int q, int n) {
	int *matrix = malloc(n * n * sizeof(int));
	int *perm = malloc(n * sizeof(int));

	for (int i = 0; i < n * n; i++)
		matrix[i] = (i % (n + 1) == 0) ? 1 : 0;
	for (int i = 0; i < n; i++)
		perm[i] = i + 1;

	for (int i = 0; i < priv_key_len; i++)
		emult(matrix, perm, priv_key[i], t_values, q, n);

	*dest_matrix = matrix;
	*dest_perm = perm;
}

void generate_signature(int **dest, int *dest_len, char *hash, int hash_len,
		int *priv_key, int priv_key_len, int *pub_key_matrix, int *pub_key_perm,
		int *t_values, int bits) {
	int n = 8;
	int q = 32;
	int a = 0, b = 0;
	//int *t_values = malloc(n * sizeof(int));
	int rand = 0;

	int *v, *v1, *v2, v_len, v1_len, v2_len;
	int *identity_perm = malloc(n * sizeof(int));

	int mark = 0;
	int *enc, enc_len;
	int *sig, sig_len;

	for (int i = 0; i < n; i++)
		identity_perm[i] = i + 1;


	// Choose a and b.
	syscall(SYS_getrandom, &rand, sizeof(int), 0);
	a = abs(rand) % (n - 3) + 2;
	do {
		syscall(SYS_getrandom, &rand, sizeof(int), 0);
		b = abs(rand) % (n - a - 1) + a + 1;
	} while (a == b);

	// In FUNKTIONSAUFRUF ÜBERGEBEN!
	a = 2;
	b = 4;
	/*
	// Choose T-values.
	for (int i = 0; i < n; i++) {
		if (i + 1 == a || i + 1 == b) {
			t_values[i] = 1;
		} else {
			syscall(SYS_getrandom, &rand, sizeof(int), 0);
			t_values[i] = abs(rand) % (q - 1) + 1;
		}
	}
	*/

	// Generate cloaking elements.
	gen_cloaking_elem(&v, &v_len, identity_perm, a, b, n, bits);
	gen_cloaking_elem(&v1, &v1_len, pub_key_perm, a, b, n, bits);
	gen_cloaking_elem(&v2, &v2_len, pub_key_perm, a, b, n, bits);

	// Encode hash into braid.
	encode(&enc, &enc_len, hash, hash_len, n);
	
	// Build signature.
	sig_len = v_len +  2 * v1_len + v2_len + 2 * priv_key_len + enc_len;
	sig = malloc(sig_len * sizeof(int));
	
	memcpy(sig, v2, v2_len * sizeof(int));
	mark += v2_len;

	for (int i = 0; i < v1_len; i++)
		sig[mark + i] = v1[v1_len - i - 1] * -1;
	mark += v1_len;

	for (int i = 0; i < priv_key_len; i++)
		sig[mark + i] = priv_key[priv_key_len - i - 1] * -1;
	mark += priv_key_len;

	memcpy(sig + mark, v, v_len * sizeof(int));
	mark += v_len;

	memcpy(sig + mark, enc, enc_len * sizeof(int));
	mark += enc_len;

	memcpy(sig + mark, priv_key, priv_key_len * sizeof(int));
	mark += priv_key_len;

	memcpy(sig + mark, v1, v1_len * sizeof(int));
	mark += v1_len;

	assert(mark == sig_len);

	// Convert signature into left canonical form.
	sig = left_canonical_form(&sig_len, sig, sig_len, n);
	printf("SIGNATURE LENGTH: %i (after BKL)\n", sig_len);
	print_braid(sig, sig_len);

	// Apply Dehornoy's algorithm to signature.
	dehornoy(&sig, &sig_len, sig, sig_len);
	printf("SIGNATURE LENGTH: %i (after Dehornoy)\n", sig_len);
	print_braid(sig, sig_len);

	*dest = sig;
	*dest_len = sig_len;
}

int verify_signature(int *sig, int sig_len, char *hash, int hash_len,
		int *pub_key_matrix, int *pub_key_perm, int *t_values, int q, int n) {
	int *enc, enc_len;
	int *msg_pub_key_matrix, *msg_pub_key_perm;

	// Encode hash into braid.
	encode(&enc, &enc_len, hash, hash_len, n);

	// Generate public key for encoded message.
	gen_pub_key(&msg_pub_key_matrix, &msg_pub_key_perm,
		enc, enc_len, t_values, q, n);

	printf("\n\n\n\n");
	print_braid(sig, sig_len);
	printf("\n\n\n\n");

	/*
	for (int i = 0; i < sig_len; i++) {
		emult(msg_pub_key_matrix, msg_pub_key_perm, sig[i], t_values, q, n);
	}
	*/

	printf("TEST ==============================\n");
	int *matrix = calloc(n * n, sizeof(int));
	// Matrix multiplication
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				matrix[i * n + j] ^= galois_mult(
					msg_pub_key_matrix[i * n + k], pub_key_matrix[k * n + j]);
				/*
				matrix[i * n + j] ^= galois_single_multiply(
					msg_pub_key_matrix[i * n + k], pub_key_matrix[k * n + j], log2(32));
					*/
			}
		}
	}


	for (int i = 0; i < sig_len; i++)
		emult(pub_key_matrix, pub_key_perm, sig[i], t_values, q, n);
	printf("Should be equal to:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%2i, ", pub_key_matrix[i * n + j]);
		}
		printf("\n");
	}

	printf("Multiplied matrix:\n");
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

void inverse_braid(int *dest, int *braid, int len) {
	for (int i = 0; i < len; i++)
		dest[i] = braid[len - i - 1];
}

/* Encodes a message given as a bit stream into a braid. The generators used are
 * picked at random, thus this function is non-deterministic. */
void encode(int **dest, int *dest_size, char *msg, int msg_size, int n) {
	// Reserve enough space for the unreduced worst case.
	int *braid = calloc(PURE_BRAID_SIZE(1, n) * 4 * msg_size * 2, sizeof(int));
	unsigned char digest = *msg;

	int mark = 0;            // Marker for the braid array index.
	int gens[4] = {0};       // Generators to use in the encoded message.
	int gen1, gen2,          // Generators for the digests, two low bits.
		power1, power2;      // Powers of the generators, two high bits.
	int reduced_middle = 0;  // Flag: Remove one reduced middle generator?
	int digest_complete = 0; // Flag: Get new digest from message?
	int duplicate;           // Flag: Generator already in use?

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

	///* TEST
	gens[0] = 1;
	gens[1] = 3;
	gens[2] = 5;
	gens[3] = 7;
	//*/

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

void gen_pure_braids(int **dest, int *dest_size, int min_len, int n) {
	int *braid = calloc(min_len + PURE_BRAID_SIZE(1, n), sizeof(int));
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
void gen_pure_braid(int *dest, int x, int n) {
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

/* Generates a private key in the braid group n with the specified fixed values
 * a and b and the specified security given in bits. */
void gen_priv_key2(int **dest, int *dest_size, int n, int bits) {
	int length = 256; // how do we calculate this?
	int *braid = calloc(length, sizeof(int));
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

/* Generates a cloaking elemtent in the braid group Bn with the specified delimiters a, b and c, with the specified number of bits in security. */
//void gen_cloaking_elem(int **dest, int *dest_size,
		//int a, int b, int n, int bits) {
void gen_cloaking_elem(int **dest, int *dest_size,
		int *perm, int a, int b, int n, int bits) {
	int mark = 0;
	//int gens_len = bits / (2 * log2f((n - 1) * (n - 2))) + 1;
	int gens_len = bits / (3 * log2f(n * (n - 1))) + 1;
	printf("Min Generators Length: %i\n", gens_len);

	// Allocate enough space for the worst case.
	// This needs to be refined!
	//int *braid = malloc(PURE_BRAID_SIZE(1, n - 1) * 4 * l * sizeof(int));
	int *braid = malloc((10 + 4 + 14 * 2 * gens_len) * sizeof(int));
	// What about the permutation part at the very start and end???
	
	int rand = 0;
	syscall(SYS_getrandom, &rand, sizeof(int), 0);
	int mid = abs(rand) % (n - 1) + 1;
	printf("Middle choice: %i\n", mid);

	// Get a and b of inverse permutation.
	int perm_a = 0, perm_b = 0;
	for (int i = 0; i < n; i++) {
		if (perm[i] == a)
			perm_a = i + 1;
		else if (perm[i] == b)
			perm_b = i + 1;
	}
	printf("Perm-A: %i, Perm-B: %i\n", perm_a, perm_b);

	// Build permuting part at start to satisfy mid -> a, mid + 1 -> b.
	if (mid <= perm_a) {
		for (int i = 0; i < perm_a - mid; i++, mark++)
			braid[mark] = perm_a - i - 1;
		for (int i = 0; i < perm_b - (mid + 1); i++, mark++)
			braid[mark] = perm_b - i - 1;
	} else if (mid < perm_b) {
		for (int i = 0; i < perm_b - (mid + 1); i++, mark++)
			braid[mark] = mid + 1 + i;
		for (int i = 0; i < mid - perm_a; i++, mark++)
			braid[mark] = perm_a + i;
	} else {
		for (int i = 0; i < (mid + 1) - perm_b; i++, mark++)
			braid[mark] = perm_b + i;
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
	int *pures = NULL;
	gen_pure_braids(&pures, &gens_size, gens_len, n);
	puts("Generated Pures:");
	print_braid(pures, gens_size);
	memcpy(braid + mark, pures, gens_size * sizeof(int));
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

//void emult(struct cb_pair *result, struct cb_pair m, int braid, int *t_values) {
void emult(int *matrix, int *perm, int artin, int *t_values, int q, int n) {
	int w = log2(q); // Prime potency for galois field arithmetic
	int abs_artin = abs(artin);
	int *result_matrix = calloc(n * n, sizeof(int));
	int *result_perm = calloc(n, sizeof(int));

	printf("Artin generator: %i\n", artin);

	// T-variable permutation and evaluation
	int t_val;
	if (artin < 0) {
		int t_i = perm[(abs_artin - 1 + 1) % n];
		//t_val = galois_inverse(t_values[t_i - 1], w);
		t_val = galois_inverse(t_values[t_i - 1]);
	} else {
		int t_i = perm[artin - 1];
		t_val = t_values[t_i - 1];
	}

	printf("Evaluated t-value: %i\n", t_val);

	// Colored Burau matrix construction
	int artin_matrix[n * n];
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
					matrix[i * n + k], artin_matrix[k * n + j]);
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

	memcpy(matrix, result_matrix, n * n * sizeof(int));
	memcpy(perm, result_perm, n * sizeof(int));
	free(result_matrix);
	free(result_perm);
}

/*
void cloaking_example() {
	int *cloak, cloak_size;
	gen_cloaking_elem(&cloak, &cloak_size, 3, 6, 8, 32);
	printf("Cloaking element size: %i\n", cloak_size);

	struct cb_pair cb;
	cb.n = 8;
	cb.q = 32;
	int cb_matrix[8*8] = {
		15, 30, 7, 18, 13, 20, 15, 31,
		10, 19, 13, 19, 6, 17, 11, 21,
		10, 14, 16, 19, 6, 17, 11, 21,
		17, 31, 21, 28, 15, 23, 2, 16,
		9, 16, 10, 13, 12, 7, 31, 20,
		9, 16, 10, 13, 21, 30, 31, 20,
		0, 0, 2, 4, 23, 0, 8, 24,
		0, 0, 0, 0, 0, 0, 0, 1,
	};
	cb.matrix = cb_matrix;
	int cb_permutation[8] = {2, 4, 3, 8, 1, 6, 7, 5};
	cb.permutation = cb_permutation;

	int t_values[8] = {28, 24, 1, 9, 26, 1, 18, 18};

	for (int i = 0; i < cloak_size; i++) {
		struct cb_pair result;
		result.n = 8;
		result.q = 32;
		int result_matrix[8*8] = {0};
		result.matrix = result_matrix;
		int result_permutation[8] = {0};
		result.permutation = result_permutation;

		emult(&result, cb, cloak[i], t_values);

		memcpy(cb.matrix, result.matrix, 8 * 8 * sizeof(int));
		memcpy(cb.permutation, result.permutation, 8 * sizeof(int));
	}
}
*/

/*
void example() {
	struct cb_pair m;
	m.n = 8;
	m.q = 32;

	int m_matrix[8*8] = {
		1, 0, 0, 0, 0, 0, 0, 0,
		0, 1, 0, 0, 0, 0, 0, 0,
		0, 0, 1, 0, 0, 0, 0, 0,
		0, 0, 0, 1, 0, 0, 0, 0,
		0, 0, 0, 0, 1, 0, 0, 0,
		0, 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0, 0, 0, 1
	};
	m.matrix = malloc(8 * 8 * sizeof(int));
	memcpy(m.matrix, &m_matrix, 8 * 8 * sizeof(int));

	int m_permutation[8] = {1, 2, 3, 4, 5, 6, 7, 8};
	m.permutation = malloc(8 * sizeof(int));
	memcpy(m.permutation, &m_permutation, 8 * sizeof(int));

	int braid1[] = {
		-4,3,-6,1,-4,-3,1,-4,-5,2,-1,-2,-2,-3,-1,5,5,-1,4,4,5,4,3,-4,-7,-3,6,5,
		-7,4,2,1,-7,-5,1,2,6,-4,-1,-2,-5,-1,4,-3,6,-3,5,1,-5,-2,4,-6,-7,-1,-1,3,
		-7,-4,-3,2,-5,2,5,-1,6,4,2,3,4,-3,2,3,-4,5,6,5,-4,2,-1,6,-7,-6
	};

	int braid2[] = {
		7,6,5,5,5,4,3,3,3,3,3,3,3,3,3,3,3,2,1,1,1,1,1,1,1,1,-2,-3,-4,5,5,5,-6,7,
		7,6,5,4,3,3,3,3,-4,5,-6,7,7,7,7,7,7,6,5,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
		3,3,3,3,3,3,3,3,3,-4,5,5,5,-6,7,7,7,7,7,7,7,7,7,7,6,5,4,3,3,-4,5,5,4,3,
		2,1,1,1,1,1,1,-2,3,3,2,1,1,1,1,1,1,-2,3,3,3,3,3,3,3,3,3,3,3,-4,5,5,5,5,
		5,5,5,-6,7,7,7,7,7,7,6,5,5,5,5,5,5,5,5,5,4,3,3,3,3,-4,5,5,5,-6,7,7,7,7,
		7,7,7,7,7,7,7,7,6,5,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,1,1,-2,3,-4,
		-5,-6,7,7,7,7,7,7,6,5,4,3,2,1,1,-2,3,-4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
		-6,7,7,7,7,7,7,7,7,6,5,4,3,2,1,1,-2,-3,-4,-5,-6,7,7,6,5,5,5,5,5,5,5,5,5,
		4,3,3,3,3,3,3,3,3,3,2,1,1,1,1,1,1,1,1,-2,-3,-4,5,5,5,5,5,5,5,5,5,5,4,3,
		3,3,3,3,3,3,3,-4,5,5,5,5,5,5,5,-6,7,7,7,7,7,7,7,7,6,5,4,3,2,1,1,1,1,1,1,
		1,1,1,1,1,1,-2,-3,-4,5,5,5,5,5,5,5,5,4,3,2,1,1,-2,-3,-4,5,5,5,5,4,3,2,1,
		1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-2,-3,-4,-5,-6,-7
	};

	int braid[] = {
		5, 4, 3, 3, -4, -5, 6, 5, 5, -6, 7, 6, 6, -7, 5, 4, 4, -5, 3, 2, 2, -3,
		6, 5, 5, -6, 7, 6, 5, 5, -6, -7, 7, 6, 6, -7, 7, 6, 6, -7, 4, 3, 3, -4,
		6, 5, 5, -6, 7, 6, 5, 4, 3, 2, 2, -3, -4, -5, -6, -7, 4, 3, 2, 2, -3, 6,
		5, 4, 3, 3, -4, -5, -6, 6, 5, 4, 4, -5, -6, 7, 6, 6, -7, 4, 3, 3, -4, 4,
		3, 2, 1, 1, -2, -3, -4, 3, 2, 1, 1, -2, -3, 6, 5, 5, -6, 6, 5, 5, -6, 5,
		4, 3, 2, 1, 1, -2, -3, -4, -5, 7, 6, 6, -7, 7, 6, 6, -7	
	};

	int tvals[8] = {28, 24, 1, 9, 26, 1, 18, 18};

	for (int i = 0; i < sizeof(braid) / sizeof(int); i++) {
		struct cb_pair result;
		result.n = 8;
		result.q = 32;
		int r_matrix[8*8] = {0};
		result.matrix = r_matrix;
		int r_permutation[8] = {0};
		result.permutation = r_permutation;

		emult(&result, m, braid[i], tvals);

		memcpy(m.matrix, result.matrix, 8 * 8 * sizeof(int));
		memcpy(m.permutation, result.permutation, 8 * sizeof(int));
	}

	free(m.matrix);
	free(m.permutation);
}
*/

/*
void example2() {
	struct cb_pair2 m;
	m.n = 4;
	m.q = 32;
	int m_matrix[4*4] = {
		1, 2, 4, 1, 
		0, 1, 2, 0, 
		0, 3, 1, 0, 
		0, 0, 0, 1
	};
	m.matrix = malloc(4 * 4 * sizeof(int));
	memcpy(m.matrix, &m_matrix, 4 * 4 * sizeof(int));
	int m_permutation[4] = {2, 3, 1, 4};
	m.permutation = malloc(4 * sizeof(int));
	memcpy(m.permutation, &m_permutation, 4 * sizeof(int));

	int braid[] = {
		3, 3, 2, 1, -3, -2, 1, 1
	};

	int tvals[4] = {3, 2, 3, 4};
	for (int i = 0; i < sizeof(braid) / sizeof(int); i++) {
		struct m_pair result;
		result.n = 4;
		result.q = 5;
		int r_matrix[4*4] = {0};
		result.matrix = r_matrix;
		int r_permutation[4] = {0};
		result.permutation = r_permutation;

		emult(&result, m, braid[i], tvals);

		memcpy(m.matrix, result.matrix, 4 * 4 * sizeof(int));
		memcpy(m.permutation, result.permutation, 4 * sizeof(int));
	}

	free(m.matrix);
	free(m.permutation);
}
*/

