CC = gcc
HWK = /c/cs323/Hwk1

all: hw1 clean

hw1: hw1.o
	${CC} -o $@ $^

clean:
	rm -f *.o
