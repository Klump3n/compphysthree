SHELL:=/bin/bash

BIN=monte_carlo

all:
	module load cuda; \
	nvcc -lm -O3 -arch sm_30 -o bin/$(BIN) src/test.cu src/geom_pbc.cu -Wno-deprecated-gpu-targets

clean:
	rm ./bin/$(BIN)
