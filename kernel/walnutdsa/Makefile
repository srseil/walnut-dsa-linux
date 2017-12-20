obj-m := walnutdsa.o
walnutdsa-objs := src/walnutdsa.o src/walnut.o src/bkl.o src/galois.o src/dehornoy.o

#KERNEL_SRC = /lib/modules/$(shell uname -r)/build
KERNEL_SRC = /lib/modules/4.14.0/build

# Hinweis: Ich musste PWD hier definieren, obwohl das nicht sein sollte.

all:
	make -C $(KERNEL_SRC) M=$(PWD) modules

clean:
	make -C $(KERNEL_SRC) M=$(PWD) clean

