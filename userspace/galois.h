#ifndef GALOIS_H
#define GALOIS_H

#include <stdint.h>

uint8_t galois_mult(unsigned int a, unsigned int b, unsigned int q);
uint8_t galois_inverse(unsigned int x, unsigned int q);

#endif

