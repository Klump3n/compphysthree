SHELL:=/bin/bash

BIN=monte_carlo_test

all:
	module load cuda; \
	nvcc -lm -O3 -arch sm_30 -o bin/$(BIN) src/test.cu src/geom_pbc.cu src/boltzmann.cu src/cases.cu -Wno-deprecated-gpu-targets

clean:
	rm ./bin/$(BIN)
