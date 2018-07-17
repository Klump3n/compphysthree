#ifndef METROPOLIS_H
#define METROPOLIS_H

#define NTRIAL 10

extern double metro_sweep(double delta, double *rnd ,int use_gpu);
extern int update_spin(int idx, double* rnd, double del, int ntrial);
void init_gpu(void);
void gather_phi(void);
extern __host__ __device__ double alocal3(spin* p, spin b, double la);
extern __host__ __device__ spin blocal(spin *p, int idx, int *nni, int nd, int nv, double ka, spin hf);

#endif
