CC=gcc
CFLAGS=-O2
#-std=c99 
# for valgrind/line debugging:
#-g

OBJECTS = aln

all: $(OBJECTS)

aln: aln.c
	$(CC) $(CFLAGS) aln.c chain.c -o aln -lz -lm

.PHONY: clean
clean:
	-rm $(OBJECTS)
