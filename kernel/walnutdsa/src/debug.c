#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "debug.h"

void print_braid(int8_t *braid, unsigned int length);

static void print_braid(int8_t *braid, unsigned int len) {
	int i;
	char seg[5] = "";
	char *complete = kcalloc(len * 5, sizeof(char), 0);
	for (i = 0; i < len; i++) {
		sprintf(seg, "%i, ", braid[i]);
		strcat(complete, seg);
	}
	printk(KERN_ALERT "%s\n", complete);
	kfree(complete);
}

