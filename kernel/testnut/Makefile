obj-m := testnut.o
testnut-objs := src/testnut.o

#KERNEL_SRC = /lib/modules/$(shell uname -r)/build
KERNEL_SRC = /lib/modules/4.14.0/build

# Hinweis: Ich musste PWD hier definieren, obwohl das nicht sein sollte.

all:
	make -C $(KERNEL_SRC) M=$(PWD) modules

clean:
	make -C $(KERNEL_SRC) M=$(PWD) clean

run:
	sudo rmmod testnut; sudo insmod ./testnut.ko

