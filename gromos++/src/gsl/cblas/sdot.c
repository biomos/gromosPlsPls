#include "../header/gsl_math.h"
#include "gsl_cblas.h"
#include "cblas.h"

float
cblas_sdot (const int N, const float *X, const int incX, const float *Y,
	    const int incY)
{
#define INIT_VAL  0.0
#define ACC_TYPE  float
#define BASE float
#include "source_dot_r.h"
#undef ACC_TYPE
#undef BASE
#undef INIT_VAL
}
