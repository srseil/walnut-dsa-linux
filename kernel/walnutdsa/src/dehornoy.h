#ifndef DEHORNOY_H
#define DEHORNOY_H

#include <linux/types.h>

void dehornoy(int8_t **dest, unsigned int *dest_len,
		int8_t *braid, unsigned int len);
void free_reduce(unsigned int *new_len, int8_t *braid, unsigned int len);

#endif

