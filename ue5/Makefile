SHELL:=/bin/bash

all:
	module load cuda; \
	nvcc -lm -O3 -arch sm_30 -o bin/run-cg src/cg.cu src/geometry.cu src/linalg.cu src/run-cg.cu -Wno-deprecated-gpu-targets

clean:
	rm ./bin/run-cg
