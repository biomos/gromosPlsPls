#include "../header/gsl_math.h"
#include "gsl_cblas.h"
#include "cblas.h"

void
cblas_zher (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	    const int N, const double alpha, const void *X, const int incX,
	    void *A, const int lda)
{
#define BASE double
#include "source_her.h"
#undef BASE
}
