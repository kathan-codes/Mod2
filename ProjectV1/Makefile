# ==============================================================================
# Makefile for Electrostatic Simulation Suite (IC103 Integrated Project)
# ==============================================================================

CC = gcc
CFLAGS = -Wall -Wextra -O2 -DINTEGRATED_BUILD
LDFLAGS = -lraylib -lopengl32 -lgdi32 -lwinmm -lm

ifeq ($(OS),Windows_NT)
    SHELL := cmd.exe
    TARGET := electric_field_simulator.exe
else
    TARGET := electric_field_simulator
endif

# Sources and Objects
SRCS = main.c \
       charge_config.c \
       electric_field.c \
       equipotential.c \
       electric_field_mapper.c \
       results_manager.c

OBJS = $(SRCS:.c=.o)

HEADERS = charge_config.h \
          electric_field.h \
          equipotential.h \
          electric_field_mapper.h \
          results_manager.h

# Default Target: Build executable
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(TARGET)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Run Target: Execute compiled application
run: $(TARGET)
	./$(TARGET)

# Clean Target: Removes generated object files and executable
clean:
ifeq ($(OS),Windows_NT)
	-@for %%F in ($(OBJS) $(TARGET)) do if exist %%F del /q /f %%F
else
	-rm -f $(OBJS) $(TARGET)
endif

.PHONY: all run clean
