#include <cuComplex.h>

#define DBL_EPSILON 2.2204460492503131e-16

#ifdef DEFINE_GLOBAL
#     define EXTERN
#else
#     define EXTERN extern
#endif

EXTERN int        ndim, nvol, *lsize, **nn;

EXTERN cuDoubleComplex *phi;

#undef EXTERN
