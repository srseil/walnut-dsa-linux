#ifndef BKL_H
#define BKL_H

#include <linux/types.h>

void left_canonical_form(int8_t **dest, unsigned int *dest_len,
	int8_t *braid, unsigned int len, unsigned int n);

#endif

