SHELL:=/bin/bash
 
all:
	module load cuda; \
	nvcc -30 -arch sm_03 bin/reduction reduction.cui -Wno-deprecated-gpu-targets

