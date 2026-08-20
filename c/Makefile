CC      ?= gcc
CFLAGS  ?= -O2 -std=gnu11 -Wall -Wextra -Wno-unused-function
CPPFLAGS += -Iinclude
LDFLAGS ?=
LDLIBS  ?= -lm

LIB_SRC := \
	src/ewcr_math.c \
	src/ewcr_signals.c \
	src/asbc.c \
	src/dwt_hybrid.c \
	src/cs_omp.c \
	src/mmc.c \
	src/svdcs.c

LIB_OBJ := $(patsubst src/%.c,build/%.o,$(LIB_SRC))
BIN_SYN := build/verify_five_codecs
BIN_FLD := build/verify_field_data

.PHONY: all test test-data clean

all: $(BIN_SYN) $(BIN_FLD)

build/%.o: src/%.c include/ewcr.h
	@mkdir -p build results
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

$(BIN_SYN): $(LIB_OBJ) build/verify_five_codecs.o
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

$(BIN_FLD): $(LIB_OBJ) build/verify_field_data.o
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

test: $(BIN_SYN)
	@mkdir -p results
	./$(BIN_SYN)

test-data: $(BIN_FLD)
	@mkdir -p results
	./$(BIN_FLD)

clean:
	rm -rf build
