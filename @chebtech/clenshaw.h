#pragma once

#ifdef HAVE_WINDOWS
#   include <intrin.h>
#else
#   include <x86intrin.h>
#endif

#include "mex.h"

// TODO: add 512-bit vectors?

#ifdef HAVE_AVX
    typedef __m256d vdouble;
#   define STRIDE 4 /* number of doubles */

#   define load_vector(_p) _mm256_load_pd(&(_p))
#   define store_vector(_p, _v) _mm256_store_pd(&(_p), (_v))
#   define broadcast_val(_s) _mm256_broadcast_sd(&(_s))

#   define mul_pd _mm256_mul_pd
#   define add_pd _mm256_add_pd
#   define sub_pd _mm256_sub_pd

#elif defined HAVE_SSE2
    typedef __m128d vdouble;
#   define STRIDE 2 /* number of doubles */

#   define load_vector(_p) _mm_load_pd(&(_p))
#   define store_vector(_p, _v) _mm_store_pd(&(_p), (_v))
#   define broadcast_val(_s) _mm_set_pd1((_s))

#   define mul_pd _mm_mul_pd
#   define add_pd _mm_add_pd
#   define sub_pd _mm_sub_pd

#else
    typedef double vdouble;
    define STRIDE 1

#   define load_vector(_p) (_p)
#   define store_vector(_p, _v) (_p) = (_v)
#   define broadcast_val(_s) (_s)

#   define mul_pd(_a, _b) ((_a) * (_b))
#   define add_pd(_a, _b) ((_a) + (_b))
#   define sub_pd(_a, _b) ((_a) - (_b))

#endif

#ifdef HAVE_FMA
#   ifdef HAVE_AVX
#       define fmadd_pd _mm256_fmadd_pd
#       define fmsub_pd _mm256_fmsub_pd
#   else
#       error "Can't compile with FMA but without AVX"
#   endif
#else
     // supports any size SIMD vectors, or scalar code.
#    define fmadd_pd(_a, _b, _c) (add_pd(mul_pd((_a), (_b)), (_c)))
#    define fmsub_pd(_a, _b, _c) (sub_pd(mul_pd((_a), (_b)), (_c)))
#endif

#define ALIGNMENT (sizeof(double) * STRIDE)
#define ALIGN_UP(_s, _a) ((size_t)(((_s) + (_a - 1)) & ~(_a - 1)))
#define ALIGN_VECTOR(_s) ALIGN_UP((_s) * sizeof(double), ALIGNMENT)
