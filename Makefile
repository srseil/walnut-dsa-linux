
main:
	gcc -std=c99 -lm -g -pg test.c

run:
	./a.out

debug:
	./a.out

valgrind:
	valgrind ./a.out

