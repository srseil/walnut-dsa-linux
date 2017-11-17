#ifndef DEHORNOY_H
#define DEHORNOY_H

#include <stdint.h>

void dehornoy(int8_t **dest, int *dest_len, int8_t *braid, int len);
void free_reduce(int *new_len, int8_t *braid, int len);

#endif

