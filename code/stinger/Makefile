CC=gcc -std=gnu9x
CXX=g++
CFLAGS=-Iinclude -Iinclude/gtri-1 -Iinclude/gtri-2 -Iinclude/gtri-3 -Iinclude/gtri-6 `pkg-config --cflags glib-2.0` -Iinclude/bc -Isrc/bc -fopenmp -g -O2
LD_LIBS= -Llib -lm -lrt `pkg-config --libs glib-2.0`

GTRI1_SRC=src/gtri-1/stinger.c src/gtri-1/stinger-deprecated.c src/gtri-1/stinger-iterator.c src/gtri-1/stinger-physmap.c src/gtri-1/stinger-utils.c src/gtri-1/timer.c src/gtri-1/x86-full-empty.c src/gtri-1/xmalloc.c src/gtri-1/stinger-physmap-new.c src/gtri-1/csv.c src/gtri-1/stinger-shared.c
GTRI1_OBJ=$(subst .c,.o,$(subst src,obj,$(GTRI1_SRC)))

BC_SRC= src/bc/bcAlgorithms.c src/bc/bcTreeDS.c src/bc/csr.c src/bc/dsUtils.c src/bc/list.c src/bc/timing_util.c src/bc/streaming_utils.c
BC_OBJ=$(subst .c,.o,$(subst src,obj,$(BC_SRC)))

GTRI2_SRC= src/gtri-2/clustering.c src/gtri-2/bfs.c
GTRI2_OBJ=$(subst .c,.o,$(subst src,obj,$(GTRI2_SRC)))

GTRI3_SRC= src/gtri-3/seed_set_exp.c
GTRI3_OBJ=$(subst .c,.o,$(subst src,obj,$(GTRI3_SRC)))

GTRI6_SRC=src/gtri-6/evicomb.c src/gtri-6/alg_convert.c
GTRI6_OBJ=$(subst .c,.o,$(subst src,obj,$(GTRI6_SRC)))

.PHONY: all
all: main

main: main.c $(GTRI1_OBJ) $(BC_OBJ) $(GTRI2_OBJ) $(GTRI3_OBJ) $(GTRI6_OBJ) 
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LD_LIBS)

sandy: sandy.c $(GTRI1_OBJ) $(BC_OBJ) $(GTRI2_OBJ) $(GTRI3_OBJ) $(GTRI6_OBJ) 
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LD_LIBS)

obj/gtri-1/x86-full-empty.o: src/gtri-1/x86-full-empty.c
	$(subst O1,O0,$(subst O2,O0,$(subst O3,O0,$(CC) $(CPPFLAGS) $(CFLAGS) -c $^ -o $@ $(LDFLAGS) $(LD_LIBS))))

obj/%.o: src/%.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS) 

.PHONY: clean
clean:
	rm -rf main obj/*.o obj/gtri-1/*.o obj/gtri-2/*.o obj/gtri-3/*.o obj/gtri-6/*.o obj/bc/*
