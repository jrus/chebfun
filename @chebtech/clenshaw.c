/*=================================================================
 * clenshaw.c
 *
 * Evaluate a polynomial in Chebyshev basis over the interval
 * [-1, 1] at a vector of evaluation points. Called from matlab
 * with two inputs X and C, representing the vector of points and
 * a matrix of coefficients, with one polynomial per column.
 * Returns a vector Y with the same number of rows as X and the
 * same number of columns as C.
 *
 * Input:   vector X, matrix C
 * Output:  matrix Y
 *
 *=================================================================*/
#include "clenshaw.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Allocate an aligned block of memory using whatever API is known
 * to the compiler. */
static inline double *
alloc_tmp(mwSize size)
{
    double * ptr = NULL;
#ifdef HAVE_POSIX_ALLOC
    posix_memalign((void **)&ptr, ALIGNMENT, ALIGN_VECTOR(size));
#elif defined HAVE_C11_ALLOC
    ptr = aligned_alloc(ALIGNMENT, ALIGN_VECTOR(size));
#elif defined HAVE_WINDOWS_ALLOC
    ptr = _aligned_malloc(ALIGN_VECTOR(size), ALIGNMENT);
#else
#error "Your C compiler can't deal with aligned allocation. Sorry!"
#endif
    return ptr;
}

/**
 * @param[out] y: output points, memory should be preallocated
 * @param[in]  x: input points
 * @param[in]  c: coefficients
 *
 * Both px and py should be pointers to SIMD-vector-size aligned arrays of
 * doubles, of length xlen. The array of coefficients pc need not be aligned.
 * Use function alloc_tmp to ensure correct alignment.
 *
 * @return:
 */
void
clenshaw(double *py,
         size_t xlen, double *px,
         size_t clen, double *pc)
{
    int i, k;
    unsigned d;
    vdouble mx, mx2, mt;
    vdouble b1, b2;
    double b, b0 = 0;
    const double two = 2.0;

    /* Degree is one less than the number of coefficients. */
    d = clen - 1;
    if (d & 1) {
        /* For odd degree polynomial, do one assignment outside the loop. */
        b = pc[d];
        d--;
    }

    /* Loop over the evaluation points X, each time taking as many points as
     * fit in a SIMD vector. Call our vector of X points mx. */
    for (i = 0; i < xlen; i += STRIDE) {
        mx = load_vector(px[i]);
        mx2 = mul_pd(mx, broadcast_val(two));
        b2 = broadcast_val(b0);
        b1 = broadcast_val(b);

        for (k = d; k > 0; k -= 2) {
            mt = fmsub_pd(mx2, b1, b2);
            b2 = add_pd(mt, broadcast_val(pc[k]));
            mt = fmsub_pd(mx2, b2, b1);
            b1 = add_pd(mt, broadcast_val(pc[k-1]));
        }

        mt = fmsub_pd(mx, b1, b2);
        b2 = add_pd(mt, broadcast_val(pc[0]));

        store_vector(py[i], b2);
    }
}

void
clenshaw_complex(double *pyr, double *pyi,
                 size_t xlen, double *pxr, double *pxi,
                 size_t clen, double *pcr, double *pci)
{
    int i, k;
    unsigned d;
    vdouble mxr, mxi, mx2r, mx2i;
    vdouble b1r, b1i, b2r, b2i;

    double br, bi, b0 = 0;
    const double two = 2.0;

    d = clen - 1;
    if (d & 1) {
        br = pcr[d];
        bi = pci[d];
        d--;
    }

    for (i = 0; i < xlen; i += STRIDE) {
        mxr = load_vector(pxr[i]);
        mxi = load_vector(pxi[i]);
        mx2r = mul_pd(mxr, broadcast_val(two));
        mx2i = mul_pd(mxi, broadcast_val(two));
        b2r = broadcast_val(b0);
        b2i = broadcast_val(b0);
        b1r = broadcast_val(br);
        b1i = broadcast_val(bi);

        /* To compute B2 = 2*X*B1 - B2 + Ck using complex arithmetic, we have:
         *
         * real(B2) = real(2*X) * real(B1) - imag(2*X) * imag(B1) - real(B2) + real(Ck)
         * imag(B2) = real(2*X) * imag(B1) + imag(2*X) * real(B1) - imag(B2) + imag(Ck)
         * etc.
         */
        for (k = d; k > 0; k -= 2) {
            b2r = sub_pd(fmadd_pd(mx2r, b1r, broadcast_val(pcr[k])),
                         fmadd_pd(mx2i, b1i, b2r));
            b2i = add_pd(fmadd_pd(mx2r, b1i, broadcast_val(pci[k])),
                         fmsub_pd(mx2i, b1r, b2i));

            b1r = sub_pd(fmadd_pd(mx2r, b2r, broadcast_val(pcr[k-1])),
                         fmadd_pd(mx2i, b2i, b1r));
            b1i = add_pd(fmadd_pd(mx2r, b2i, broadcast_val(pci[k-1])),
                         fmsub_pd(mx2i, b2r, b1i));
        }

        b2r = sub_pd(fmadd_pd(mxr, b1r, broadcast_val(pcr[0])),
                     fmadd_pd(mxi, b1i, b2r));
        b2i = add_pd(fmadd_pd(mxr, b1i, broadcast_val(pci[0])),
                     fmsub_pd(mxi, b1r, b2i));

        store_vector(pyr[i], b2r);
        store_vector(pyi[i], b2i);
    }
}

/* Call the matlab function full to convert a sparse -> full array. */
mxArray *
sparseToFull(const mxArray *in) {
    mxArray *plhs[1];
    const mxArray *prhs[1] = {in};
    mexCallMATLAB(1, plhs, 1, (mxArray **)prhs, "full");
    return plhs[0];
}

/* Main function. */
void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *px, *pc;
    mxArray *py;
    mwSize x_nrow, c_nrow, c_ncol, i;
    mxComplexity x_complexity, c_complexity, y_complexity;
    double *tmpx = NULL, *tmpxi = NULL, *tmpy = NULL, *tmpyi = NULL;
    double *colyr, *colyi, *colcr, *colci, *tmpci = NULL;

    /* Check for proper number of arguments. */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("CHEBFUN:CHEBTECH:clenshaw:inputs",
            "Clenshaw requires two inputs, a vector of evaluation points X "
            "and a matrix of coefficients C.");
    }
    /* Make sure there is no more than one output. */
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("CHEBFUN:CHEBTECH:clenshaw:outputs",
            "Too many output arguments. Clenshaw returns only one output.");
    }
    /* Give the input pointers more descriptive names. */
    px = prhs[0];
    pc = prhs[1];

    /* Make sure that X and C are both matrices of doubles. */
    if (!mxIsDouble(px) || !mxIsDouble(pc)) {
        mexErrMsgIdAndTxt("CHEBFUN:CHEBTECH:clenshaw:inputtype",
            "Inputs to Clenshaw must be of type double.");
    }
    /* Check that X is a column vector, and that C is a matrix. */
    if (mxGetN(px) != 1 || mxGetNumberOfDimensions(px) > 2) {
        mexErrMsgIdAndTxt("CHEBFUN:CHEBTECH:clenshaw:xDim",
            "Evaluation points X must be a column vector.");
    }
    if (mxGetNumberOfDimensions(pc) > 2) {
        mexErrMsgIdAndTxt("CHEBFUN:CHEBTECH:clenshaw:cDim",
            "Coefficients C must be either a column vector or a "
            "two-dimensional matrix.");
    }

    x_nrow = mxGetM(px);
    c_nrow = mxGetM(pc);
    c_ncol = mxGetN(pc);

    if (c_ncol > 0 && c_nrow == 0) {
        mexErrMsgIdAndTxt("CHEBFUN:CHEBTECH:clenshaw:cDim",
            "Coefficients C cannot be an empty 0-by-n matrix with n > 0.");
    }
    /* If X or C is empty, return a correctly shaped empty array. */
    if (x_nrow == 0 || c_ncol == 0) {
        plhs[0] = mxCreateDoubleMatrix(x_nrow, c_ncol, mxREAL);
        return;
    }

    /* If either X or C is a sparse matrix, call out to Matlab to convert it
     * to a full matrix. */
    if (mxIsSparse(px)) {
        px = sparseToFull(px);
    }
    if (mxIsSparse(pc)) {
        pc = sparseToFull(pc);
    }

    /* Test for complex inputs. */
    x_complexity = mxIsComplex(px) ? mxCOMPLEX : mxREAL;
    c_complexity = mxIsComplex(pc) ? mxCOMPLEX : mxREAL;
    y_complexity = (x_complexity == mxCOMPLEX || c_complexity == mxCOMPLEX) ?
        mxCOMPLEX : mxREAL;

    /* Allocate output matrix, and put its pointer on the output array. */
    py = mxCreateDoubleMatrix(x_nrow, c_ncol, y_complexity);
    plhs[0] = py;

    /* Create two temporary arrays which are properly memory aligned,
     * to be used for the input and output to the clenshaw function. */
    tmpx = alloc_tmp(x_nrow);
    if (tmpx == NULL) { goto oom_error; }

    tmpy = alloc_tmp(x_nrow);
    if (tmpy == NULL) { goto oom_error; }

    /* Get a pointer to the beginning of the arrays for C and I.
     * Below, we will loop through the columns, updating these to point
     * at the beginning of each column in turn. */
    colcr = mxGetPr(pc);
    colyr = mxGetPr(py);
    if (y_complexity == mxCOMPLEX) {
        colyi = mxGetPi(py);

        if (c_complexity == mxCOMPLEX) {
            colci = mxGetPi(pc);
        /* If X is complex and C is real, allocate an array of zeros as the
         * imaginary part of C */
        } else {
            tmpci = calloc(mxGetNumberOfElements(pc), sizeof(double));
            if (tmpci == NULL) { goto oom_error; }
            colci = tmpci;
        }
    }

    /* If X is real, then we can use a real-valued version of Clenshaw. */
    if (x_complexity == mxREAL) {
        /* Copy the values of X into our memory aligned array. */
        memcpy(tmpx, mxGetPr(px), x_nrow * sizeof(double));

        /* Handle real part. For each column of coefficients, run clenshaw,
         * then copy the result values to the proper column of Y. */
        for (i = 0; i < c_ncol; ++i, colyr += x_nrow, colcr += c_nrow) {
            clenshaw(tmpy, x_nrow, tmpx, c_nrow, colcr);
            memcpy(colyr, tmpy, x_nrow * sizeof(double));
        }
        /* If C is complex, separately handle imaginary part. */
        if (c_complexity == mxCOMPLEX) {
            for (i = 0; i < c_ncol; ++i, colyi += x_nrow, colci += c_nrow) {
                clenshaw(tmpy, x_nrow, tmpx, c_nrow, colci);
                memcpy(colyi, tmpy, x_nrow * sizeof(double));
            }
        }
    }

    /* If X is complex, then we need to call a complex-valued version
     * of clenshaw. */
    else {
        /* First allocate two more temporary arrays, for storing the
         * imaginary parts of our inputs and outputs. */
        tmpxi = alloc_tmp(x_nrow);
        if (tmpxi == NULL) { goto oom_error; }
        tmpyi = alloc_tmp(x_nrow);
        if (tmpyi == NULL) { goto oom_error; }

        /* Copy the values of X into our memory-aligned arrays. */
        memcpy(tmpx, mxGetPr(px), x_nrow * sizeof(double));
        memcpy(tmpxi, mxGetPi(px), x_nrow * sizeof(double));

        /* For each column of coefficients, run clenshaw_complex,
         * then copy the result values to the proper column of Y. */
        for (i = 0; i < c_ncol; ++i,
                colyr += x_nrow, colyi += x_nrow,
                colcr += c_nrow, colci += c_nrow) {
            clenshaw_complex(tmpy, tmpyi,
                             x_nrow, tmpx, tmpxi,
                             c_nrow, colcr, colci);
            memcpy(colyr, tmpy, x_nrow * sizeof(double));
            memcpy(colyi, tmpyi, x_nrow * sizeof(double));
        }

        /* Clean up temporary arrays. */
        free(tmpxi);
        free(tmpyi);
        free(tmpci);
    }

    /* Clean up temporary arrays. */
    free(tmpx);
    free(tmpy);

    return;

/* If we run out of memory, clean up and then throw an error. */
oom_error:
    free(tmpx);
    free(tmpy);
    free(tmpxi);
    free(tmpyi);
    free(tmpci);
    /* Note, mexErrMsgIdAndTxt automatically frees any memory allocated by
     * calls to mxCreateDoubleMatrix. */
    mexErrMsgIdAndTxt("CHEBFUN:CHEBTECH:clenshaw:outofmemory",
        "Insufficient Memory.");
}
