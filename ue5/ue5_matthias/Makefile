SHELL:=/bin/bash

all:
	module load cuda; \
	nvcc -lm -O3 -arch sm_30 -o run-cg cg.cu geometry.cu linalg.cu run-cg.cu -Wno-deprecated-gpu-targets

clean:
	rm ./bin/reduction
	rm ./measurements/*.txt
	rm ./src/tmpfile.txt
