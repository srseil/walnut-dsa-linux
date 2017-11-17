#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "debug.h"

void print_braid(int8_t *braid, int length);
void print_permutation(int8_t *braid, int length);

void print_braid(int8_t *braid, int length) {
	int i;
	for (i = 0; i < length; i++)
		printf("%i ", braid[i]);
	printf("\n");
}

void print_permutation(int8_t *braid, int length) {
	int perm[] = {1, 2, 3, 4, 5, 6, 7, 8};
	int i;

	for (i = 0; i < length; i++) {
		int k = abs(braid[i]);
		int temp = perm[k - 1];
		perm[k - 1] = perm[k];
		perm[k] = temp;
	}

	for (i = 0; i < 8; i++)
		printf("%i ", perm[i]);
	printf("\n");
}

