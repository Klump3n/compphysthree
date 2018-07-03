SHELL:=/bin/bash

BIN=nice_binary

all:
	module load cuda; \
	nvcc -lm -O3 -arch sm_30 -o bin/$(BIN) src/run.cu src/geom_pbc.cu src/randgpu.cu src/metropolis.cu src/spin.cu -Wno-deprecated-gpu-targets

clean:
	rm ./bin/$(BIN)
