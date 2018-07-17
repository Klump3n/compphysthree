#include "spin.h"

#define DBL_EPSILON 2.2204460492503131e-16

#ifdef DEFINE_GLOBAL
#     define EXTERN
#else
#     define EXTERN extern
#endif

EXTERN int ndim, nvol, *lsize, **nn, **eo;
EXTERN double lambda, kappa;
EXTERN spin h;

EXTERN spin *phi;

#undef EXTERN
