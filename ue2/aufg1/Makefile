cg_nvcc: cg.c
	nvcc -lm src/cg.c -o bin/conj_grad -Wno-deprecated-gpu-targets

cg_nvcc_debug: cg.c
	nvcc -lm src/cg.c -o bin/conj_grad -g -Wno-deprecated-gpu-targets

cg_gcc: cg.c
	gcc -lm src/cg.c -o bin/conj_grad

cg_gcc_debug: cg.c
	gcc -lm src/cg.c -o bin/conj_grad -g
