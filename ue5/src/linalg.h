#ifndef LINALG_H
#define LINALG_H

extern void random_vector(double *p);
extern double norm_sqr(double *v);
extern double vector_prod(double *v, double *w);
extern void assign_v2v(double *v, double *w);
extern void mul_add(double *v, double a, double *w);
extern void update_p(double *r, double b, double *p);
extern void laplace_2d(double *w, double *v);

extern __global__ void reduceUnrolling (double *g_idata, double *g_odata, unsigned int n);
extern __global__ void assign_v2v_gpu(double *v, double *w, int nx, int ny);
extern __global__ void mul_add_gpu(double *v, double a, double *w, int nx, int ny);
extern __global__ void update_p_gpu(double *r, double b, double *p, int nx, int ny);
extern __global__ void laplace_2d_gpu(double *w, double *v, int nx, int ny);


/* our stuff */
extern double vector_add(double *v, const int n);
extern void norm_sqr_gpu(double *res, double *d_r, double *d_intmed1, double *d_intmed2, int Nx, int Ny);
extern __global__ void vector_multiply_entries_gpu(double *res, double *v, double *w, int nx, int ny);
extern __global__ void vector_square_entries_gpu(double *v, int nx, int ny);
extern void vector_prod_gpu(double *res, double *d_v, double *d_w, double *d_intmed1, double *d_intmed2, int Nx, int Ny);

#endif
