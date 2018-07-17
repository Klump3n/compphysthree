#ifndef SPIN_H
#define SPIN_H

#include <cuComplex.h>

typedef cuDoubleComplex spin;

#define make_spin(a,b) ( \
   make_cuDoubleComplex(a,b) \
)

extern double action(void);
extern void random_cnfg(void);
extern spin magnet(void);
extern double alocal(int idx);
extern double alocal2(int idx);
extern spin compute_b(int idx);

#endif
