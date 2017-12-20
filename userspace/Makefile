main:
	gcc -Wall -std=gnu89 -lkcapi -lm -g main.c galois.c debug.c

run:
	./a.out

debug:
	gdb ./a.out

valgrind:
	valgrind --leak-check=full ./a.out

