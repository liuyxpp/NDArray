#CC=icpc -O2
CC=g++ -O2
#CC=g++ -g
#CC=g++ -pg -O2
CXXFLAGS=-I../lyx
LDFLAGS=
LIBS=-lm

all: CNDArray.o
	ar r libndarray.a $^

CNDArray.o:CNDArray.cpp CNDArray.h
	$(CC) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o libndarray.a
