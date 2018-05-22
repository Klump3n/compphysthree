SHELL:=/bin/bash
 
all:
	module load cuda; \
	nvcc -O3 -arch sm_30 -o bin/reduction src/reduction.cu -Wno-deprecated-gpu-targets

