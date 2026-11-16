FC = gfortran
FFLAGS = -O2 -Wall

SRC_DIR = fortran
BUILD_DIR = build

SOURCES = $(wildcard $(SRC_DIR)/*.f90)

PROGRAMS = $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%,$(SOURCES))

.PHONY: all clean list

all: $(PROGRAMS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -o $@ $<

list:
	@printf '%s\n' $(PROGRAMS)

clean:
	rm -rf $(BUILD_DIR)
