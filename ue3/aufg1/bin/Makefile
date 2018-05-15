SHELL:=/bin/bash

all:
	module load cuda; \
	nvcc -o addArrayTwoDim addArrayHost.cu -Wno-deprecated-gpu-targets
