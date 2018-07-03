SHELL:=/bin/bash

BIN=mc_mag

all:
	module load cuda; \
	nvcc -lm -O3 -arch sm_30 -o bin/$(BIN) src/run.cu src/geom_pbc.cu src/randgpu.cu src/metropolis.cu src/spin.cu src/stat5.cu -Wno-deprecated-gpu-targets

clean:
	rm ./bin/$(BIN)
