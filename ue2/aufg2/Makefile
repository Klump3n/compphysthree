SHELL:=/bin/bash

all:
	module load cuda; \
	nvcc src/addArrayHost.cu -o bin/addArrayTest -Wno-deprecated-gpu-targets
