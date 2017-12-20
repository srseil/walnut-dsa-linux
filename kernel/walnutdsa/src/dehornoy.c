#include <linux/slab.h>
#include <linux/string.h>

#include "dehornoy.h"

void dehornoy(int8_t **dest, unsigned int *dest_len,
		int8_t *braid, unsigned int len);
void free_reduce(unsigned int *new_len, int8_t *braid, unsigned int len);

static int sig(int x);


static int sig(int x) {
	return x > 0 ? 1 : (x < 0 ? -1 : 0);
}

void free_reduce(unsigned int *new_len, int8_t *braid, unsigned int len) {
	int i, j;

	for (i = 0; i < len - 1 && braid[i] != -1 * braid[i + 1]; i++);

	for (j = i + 2; j < len; j++) {
		if (j + 1 < len && braid[j] == -1 * braid[j + 1])
			j++;
		else if (i - 1 >= 0 && braid[i - 1] == -1 * braid[j])
			i--;
		else
			braid[i++] = braid[j];
	}

	*new_len = (i < len - 1) ? i : len;
}

void dehornoy(int8_t **dest, unsigned int *dest_len,
		int8_t *braid, unsigned int len) {
	int i, j;

	/* Free-reduce once at the beginning. */
	free_reduce(&len, braid, len);

	for (;;) {
		int8_t *handle;
		unsigned int handle_len;
		int start_sig;
		int main_gen;
		int off = 0;
		int search_start = 0, search_end = len;
		int start = -1, end = -1;
		int handle_found = 1;

		/* Find the left-most nested handle. */
		while (handle_found) {
			handle_found = 0;
			for (i = search_start; i < search_end - 2; i++) {
				for (j = i + 1; j < search_end; j++) {
					if (braid[j] == braid[i]
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
			/* No handle found. */
			break;
		}

		/* There does not seem to be an practical upper bound for the maximum
		 * length of a reduced handle, so we just pick 10 times the size from
		 * before. Maybe this would be a good place to use a dynamically re-
		 * sized array;
		 */
		handle = kcalloc(10 * (end - start + 1), sizeof(int8_t), 0);
		start_sig = sig(braid[start]);
		main_gen = abs(braid[start]);

		/* Reduce handle. */
		j = 1;
		for (i = 1; i < end - start; i++) {
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
		handle_len = j + 2 * off - 1;
		if (start + handle_len - 1 != end) {
			memmove(braid + start + handle_len, braid + end + 1,
					(len - end - 1) * sizeof(int8_t));
			len += start + handle_len - end - 1;
		}
		memcpy(braid + start, handle, handle_len * sizeof(int8_t));
		kfree(handle);

		free_reduce(&len, braid, len);
	}

	*dest = braid;
	*dest_len = len;
}

