CC := gcc
PCC := mpicc

ifeq ($(debug),true)
	CFLAGS := -g -O0
else
	CFLAGS := -O2
endif

BIN_PATH := bin
OBJ_PATH_BASE := bin/objs
OBJ_PATH_SERIAL := $(OBJ_PATH_BASE)/serial
OBJ_PATH_PAR := $(OBJ_PATH_BASE)/parallel
SRC_PATH_BASE := src
SRC_PATH_SERIAL := $(SRC_PATH_BASE)/serial
SRC_PATH_PAR := $(SRC_PATH_BASE)/parallel
LIBS=-lm
TARGET_NAME := wave_eq
TARGET_NAME_SERIAL := $(TARGET_NAME)_s
TARGET_NAME_PAR := $(TARGET_NAME)_p
TARGET_SERIAL := $(BIN_PATH)/$(TARGET_NAME_SERIAL)
TARGET_PAR := $(BIN_PATH)/$(TARGET_NAME_PAR)

rwildcard=$(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2) $(filter $(subst *,%,$2),$d))

SRC_SERIAL := $(call rwildcard,$(SRC_PATH_SERIAL),*.c)
OBJ_SERIAL := $(patsubst %.c,%.o,$(SRC_SERIAL:$(SRC_PATH_SERIAL)/%=$(OBJ_PATH_SERIAL)/%))
SRC_PAR := $(call rwildcard,$(SRC_PATH_PAR),*.c)
OBJ_PAR := $(patsubst %.c,%.o,$(SRC_PAR:$(SRC_PATH_PAR)/%=$(OBJ_PATH_PAR)/%))

CFLAGS_SERIAL = $(CFLAGS) $(foreach d,$(dir $(SRC_SERIAL)), -I$d)
CFLAGS_PAR = $(CFLAGS) $(foreach d,$(dir $(SRC_PAR)), -I$d)


default: makedir all

# non-phony targets
$(BIN_PATH)/$(TARGET_NAME_SERIAL): $(OBJ_SERIAL)
	$(CC) $(CFLAGS_SERIAL) -o $@ $(OBJ_SERIAL) $(LIBS)

$(OBJ_PATH_SERIAL)/%.o: $(SRC_PATH_SERIAL)/%.c
	$(CC) $(CFLAGS_SERIAL) -c -o $@ $<

$(BIN_PATH)/$(TARGET_NAME_PAR): $(OBJ_PAR)
	$(PCC) $(CFLAGS_PAR) -o $@ $(OBJ_PAR) $(LIBS)

$(OBJ_PATH_PAR)/%.o: $(SRC_PATH_PAR)/%.c
	$(PCC) $(CFLAGS_PAR) -c -o $@ $<

# phony rules
.PHONY: makedir
makedir:
	@mkdir -p $(OBJ_PATH_SERIAL) $(OBJ_PATH_PAR)

.PHONY: all
all: $(TARGET_SERIAL) $(TARGET_PAR)

.PHONY: clean
clean:
	@rm -fR $(TARGET_SERIAL) $(OBJ_SERIAL) $(TARGET_PAR) $(OBJ_PAR)
	@echo Project cleaned.

.PHONY: run_serial
run_serial: $(TARGET_SERIAL)
	@./$^

.PHONY: run_par
run_par: $(TARGET_PAR)
	@mpirun -n 5 ./$^
