SHELL:=/bin/bash

all:
	module load cuda
	nvcc -o addArrayTest addArrayHost.cu -Wno-deprecated-gpu-targets
