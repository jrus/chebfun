#pragma once

#include <x86intrin.h>
#include "mex.h"

#define STRIDE 4 /* number of doubles */
#define ALIGNMENT (sizeof(double) * STRIDE)

#define ALIGN_UP(_s, _a) ((size_t)(((_s) + (_a - 1)) & ~(_a - 1)))
#define ALIGN_VECTOR(_s) ALIGN_UP((_s) * sizeof(double), ALIGNMENT)

typedef __m256d vdouble;

#define load_pd _mm256_load_pd
#define store_pd _mm256_store_pd
#define broadcast_sd(_s) _mm256_broadcast_sd(&(_s))

#define add_pd _mm256_add_pd
#define sub_pd _mm256_sub_pd
#define mul_pd _mm256_mul_pd
#define fmadd_pd _mm256_fmadd_pd
#define fmsub_pd _mm256_fmsub_pd

static inline double *alloc_tmp(mwSize size);
void clenshaw(double *py,
              size_t xlen, double *px,
              size_t clen, double *pc);
void clenshaw_complex(double *pyr, double *pyi,
                      size_t xlen, double *pxr, double *pxi,
                      size_t clen, double *pcr, double *pci);

mxArray *sparseToFull(const mxArray *in);
