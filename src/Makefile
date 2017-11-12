CC=gcc
#CC=gcc-5.2.0
#CFLAGS=-Wall -O0 -c -DUSE_LCPX -g
CFLAGS=-Wall -Wno-unused-variable -O3 -c -DUSE_LCPX -g -fopenmp -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free #-DUSELEXICOGRAPHICAL #-DWITHGPERFTOOLS
INCLUDES=-I. -I/usr/local/include
LDFLAGS=-lm -L/usr/local/lib -fopenmp -ltcmalloc #-L/home/ajmedlar/local/lib #-lprofiler
SOURCES=main.c \
        index.c \
        demux.c \
        seq.c \
        fasta.c \
        utils.c \
        search.c \
        db.c \
        engine.c \
        rbtree.c \
        vector.c \
        hit.c \
        suffix.c \
        query.c \
        work_queue.c \
        priority_queue.c \
        alignment.c \
        ssw.c \
        ssw_wrapper.c \
        divsufsort.c \
        trsort.c \
        sssort.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=topaz

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

