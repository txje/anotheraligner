CC=gcc
CFLAGS=-O2
#-std=c99 

OBJECTS = aln

all: $(OBJECTS)

aln: aln.c
	$(CC) $(CFLAGS) aln.c chain.c -o aln -lz

.PHONY: clean
clean:
	-rm $(OBJECTS)
