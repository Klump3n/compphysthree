SHELL:=/bin/bash

all:
	scripts/factor_numbers.py
	mv grid_parameters scripts
	module load cuda; \
	nvcc -o bin/addArrayTwoDim src/addArrayHost.cu -Wno-deprecated-gpu-targets
