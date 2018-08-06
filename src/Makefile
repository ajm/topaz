CC=gcc
CFLAGS=-Wall -Wno-unused-variable -O3 -c -g #-DUSELEXICOGRAPHICAL
INCLUDES=-I. -I./seg
LDFLAGS=-lm -lpthread
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
        sssort.c \
        seg_wrapper.c \
        seg/seg.c \
        seg/genwin.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=topaz

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@ 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

