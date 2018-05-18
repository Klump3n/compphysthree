SHELL:=/bin/bash

all:
	module load cuda; \
	nvcc -lm -o bin/laplace src/laplaceGPU.cu -Wno-deprecated-gpu-targets
