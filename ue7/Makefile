SHELL:=/bin/bash

BIN=spin_update

all:
	module load cuda; \
	nvcc -lm -O3 -arch sm_30 -o bin/$(BIN) src/action.cu src/geom_pbc.cu src/eigener_code.cu src/randgpu.cu src/tests.cu -Wno-deprecated-gpu-targets

clean:
	rm ./bin/$(BIN)
