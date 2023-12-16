CC := gcc

ifeq ($(release),true)
	CFLAGS := -O2
else
	CFLAGS := -g -O0
endif

BIN_PATH := bin
OBJ_PATH_BASE := bin/objs
OBJ_PATH_SERIAL := $(OBJ_PATH_BASE)/serial
SRC_PATH_BASE := src
SRC_PATH_SERIAL := $(SRC_PATH_BASE)/serial
TARGET_NAME := wave_eq
TARGET_NAME_SERIAL := $(TARGET_NAME)_s
TARGET := $(BIN_PATH)/$(TARGET_NAME)
TARGET_SERIAL := $(BIN_PATH)/$(TARGET_NAME_SERIAL)

rwildcard=$(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2) $(filter $(subst *,%,$2),$d))

SRC_SERIAL := $(call rwildcard,$(SRC_PATH_SERIAL),*.c)
OBJ_SERIAL := $(patsubst %.c,%.o,$(SRC_SERIAL:$(SRC_PATH_SERIAL)/%=$(OBJ_PATH_SERIAL)/%))

CFLAGS += $(foreach d,$(dir $(SRC_SERIAL)), -I$d)


default: makedir all

# non-phony targets
$(TARGET_SERIAL): $(OBJ_SERIAL)
	$(CC) $(CFLAGS) -o $@ $(OBJ_SERIAL)

$(OBJ_PATH_SERIAL)/%.o: $(SRC_PATH_SERIAL)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

# phony rules
.PHONY: makedir
makedir:
	@mkdir -p $(OBJ_PATH_SERIAL)

.PHONY: all
all: $(TARGET_SERIAL)

.PHONY: clean
clean:
	@rm -fR $(TARGET_SERIAL) $(OBJ_PATH_SERIAL)
	@echo Project cleaned.

.PHONY: run_serial
run_serial: $(TARGET_SERIAL)
	@./$^
