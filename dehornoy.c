#include "debug.c"

void free_reduce(int *new_len, int *braid, int len);

void test_shit() {
	int braid[] = {
		1, -1, 5, 5, -5, 2, -3, 2, 3, -3, 4, -5, 5, 6, -6, 7
	};
	int len = 0;
	free_reduce(&len, braid, sizeof(braid) / sizeof(int));

	for (int i = 0; i < len; i++)
		printf("%i ", braid[i]);
	printf("\n");
}

int sig(int x) {
	return x > 0 ? 1 : (x < 0 ? -1 : 0);
}

void free_reduce(int *new_len, int *braid, int len) {
	int i = 0;
	for (; i < len - 1 && braid[i] != -1 * braid[i + 1]; i++);

	int mark = 0;
	int arr[1000] = {0};
	for (int j = i + 2; j < len; j++) {
		if (j + 1 < len && braid[j] == -1 * braid[j + 1]) {
			arr[mark++] = j;
			j++;
		} else
			braid[i++] = braid[j];
	}
	
	printf("FRlen: %i, FRi: %i\n", len, i);
	if ((len - i) % 2 != 0)
		print_braid(arr, mark);
	*new_len = (i < len - 1) ? i : len;
}

void dehornoy(int **dest, int *dest_len, int *braid, int len) {
	// Make sure that this braid has enough space for the reduction expansions!
#if 0
	int braid_inline[] = {
		//1, 0, 0, 2, 0, 0, 0, -2, 3, 3, 0, -3, 4, 5, -5
		//
		-1, -2, 1, 3, -2, -3, -2, 1, -3, 2, 1, 1
		//0, -1, 0, -2, 0, 0, 1, 3, -2, 0, -3, -2, 1, -3, 2, 1, 1
		// -> 2, -1, -2, 3, -2, -3, -2, 1, -3, 2, 1, 1

		//2, -1, -2, 3, -2, -3, -2, 1, -3, 2, 1, 1
		//2, 0, -1, -2, 0, 0, 3, -2, 0, -3, -2, 1, 0, -3, 2, 1, 1
		
		//2, 2, -1, 3, 2, -1, -2, -3, -1, 3, -2, -3, 1, 1
		//2, 2, 0, 0, -1, 3, 0, 2, -1, -2, 0, -3, -1, 0, 3, 0, -2, -3, 1, 1

		//2, 2, -1, 3, 2, -1, -2, -3, -1, -2, -3, 2, 1, 1

		//2, 2, -1, 3, 2, -1, 2, -1, -2, -3 // no handle
		//
		//ENDE:
		// 2, 2, -1, 3, 2, -1, 2, -1, -2, -3
	};
	len = sizeof(braid_inline) / sizeof(int);
	braid = calloc(2 * len, sizeof(int));
	memcpy(braid, braid_inline, sizeof(braid_inline));
	//braid = braid_inline;
#endif


	int orig = len;
	int cnt = 0;



	/*
	int *temp = braid;
	braid = calloc(10 * len, sizeof(int));
	memcpy(braid, temp, len * sizeof(int));
	*/


	// Free-reduce once at the beginning.
	free_reduce(&len, braid, len);

	while (true) {
		// Find the left-most nested handle.
		int search_start = 0, search_end = len;
		int start = -1, end = -1;
		int handle_found = 1;
		while (handle_found) {
			handle_found = 0;
			for (int i = search_start; i < search_end - 2; i++) {
				for (int j = i + 1; j < search_end; j++) {
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
			printf("No handle found.\n");
			break;
		}

		// 2 * too short for handle?
		int *handle = calloc(10 * (end - start + 1), sizeof(int));

		int start_sig = sig(braid[start]);
		int main_gen = abs(braid[start]);
		//printf("start_sig: %i, main_gen: %i\n", start_sig, main_gen);
		int off = 0;

		// Reduce handle.
		int j = 1;
		for (int i = 1; i < end - start; i++) {
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


		printf("start: %i, end: %i, len: %i", start, end, end - start + 1);
		printf("Pre-reduced handle: ");
		print_braid(braid + start, end - start + 1);
		int handle_len = j + 2 * off - 1;
		printf("handle_len: %i\n", j + 2 * off - 1);
		printf("Handle: ");
		print_braid(handle, handle_len);
		printf("LENGTH: %i\n", len);

		/*
		for (int i = start; i <= end; i++)
			braid[i] = 0;
			*/
		//printf("Before: ");
		//print_braid(braid, len);

		if (start + handle_len - 1 != end) {
		//if (start + handle_len > end) {
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
		//for (int i = 0; i < 10; i++)


		{
			int per[] = {1, 2, 3, 4, 5, 6, 7, 8};
			for (int i = 0; i < len; i++) {
				int k = abs(braid[i]);
				int temp = per[k - 1];
				per[k - 1] = per[k];
				per[k] = temp;
			}

			printf("PERMUTATION: ");
			for (int i = 0; i < 8; i++)
				printf("%i ", per[i]);
			printf("\n");
		}

		free_reduce(&len, braid, len);

		printf("Relevant section: ");
		print_braid(braid + start - 1, handle_len + 2);

		/*
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
		*/

		/*
		   if (++cnt == 20) {
		   printf("\nBEFORE:\n");
		   print_braid(braid, len);
		   }

		   5,6
		   2,4
		   0,2
		   */
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

		{
			int per[] = {1, 2, 3, 4, 5, 6, 7, 8};
			for (int i = 0; i < len; i++) {
				int k = abs(braid[i]);
				int temp = per[k - 1];
				per[k - 1] = per[k];
				per[k] = temp;
			}

			printf("PERMUTATION: ");
			for (int i = 0; i < 8; i++)
				printf("%i ", per[i]);
			printf("\n");
		}
	
	}

		/*
	int *result = malloc(len * sizeof(int));
	int mark = 0;
	for (int i = 0; i < len; i++) {
		if (braid[i] != 0)
			result[mark++] = braid[i];
	}
	print_braid(braid, len);
	*dest = result;
	*dest_len = mark;
	*/

	print_braid(braid, len);
	*dest = braid;
	*dest_len = len;
}

