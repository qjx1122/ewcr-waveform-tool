CC      ?= gcc
CFLAGS  ?= -O2 -std=gnu11 -Wall -Wextra -Wno-unused-function
CPPFLAGS += -Iinclude
LDFLAGS ?=
LDLIBS  ?= -lm

SRC := \
	src/ewcr_math.c \
	src/ewcr_signals.c \
	src/asbc.c \
	src/dwt_hybrid.c \
	src/cs_omp.c \
	src/mmc.c \
	src/svdcs.c \
	src/verify_five_codecs.c

OBJ := $(patsubst src/%.c,build/%.o,$(SRC))
BIN := build/verify_five_codecs

.PHONY: all test clean

all: $(BIN)

build/%.o: src/%.c include/ewcr.h
	@mkdir -p build results
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) $(LDLIBS) -o $@

test: $(BIN)
	@mkdir -p results
	./$(BIN)

clean:
	rm -rf build
