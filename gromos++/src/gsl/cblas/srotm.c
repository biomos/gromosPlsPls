#include "../header/gsl_math.h"
#include "gsl_cblas.h"
#include "cblas.h"

void
cblas_srotm (const int N, float *X, const int incX, float *Y, const int incY,
	     const float *P)
{
#define BASE float
#include "source_rotm.h"
#undef BASE
}
