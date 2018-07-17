#ifndef RANDOM_H
#define RANDOM_H
#include <curand_kernel.h> // CURAND Bibliothek header-Datei

extern double* randgpu(unsigned int N);
extern double* devRandgpu(unsigned int N);
extern curandState* get_curand_states(unsigned int N);

#endif
