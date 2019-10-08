CC=gcc
CFLAGS=-O2
#-std=c99 
# for valgrind/line debugging:
#-g

ifeq ($(sse2),)
	CFLAGS += -march=native
endif
ifeq ($(avx2),)
	CFLAGS += -mavx2
endif

OBJECTS = aln

all: $(OBJECTS)

aln: aln.c
	$(CC) $(CFLAGS) aln.c chain.c ksw2_extd2_sse.c -o aln -lz -lm

.PHONY: clean
clean:
	-rm $(OBJECTS)
