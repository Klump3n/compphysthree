SHELL:=/bin/bash

all:
	module load cuda; \
	nvcc -lm -o laplace laplaceGPU.cu -Wno-deprecated-gpu-targets
