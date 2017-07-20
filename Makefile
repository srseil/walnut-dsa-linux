
main:
	gcc -std=c99 -lm -g test.c galois.c galois.h

run:
	./a.out

debug:
	./a.out

valgrind:
	valgrind ./a.out

