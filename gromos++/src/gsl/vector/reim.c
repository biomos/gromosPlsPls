#include "../../../config.h"
#include <stdlib.h>
#include "gsl_vector.h"

#include "view.h"

#define BASE_GSL_COMPLEX_LONG
#include "../header/templates_on.h"
#include "reim_source.c"
#include "../header/templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../header/templates_on.h"
#include "reim_source.c"
#include "../header/templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../header/templates_on.h"
#include "reim_source.c"
#include "../header/templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define USE_QUALIFIER
#define QUALIFIER const

#define BASE_GSL_COMPLEX_LONG
#include "../header/templates_on.h"
#include "reim_source.c"
#include "../header/templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../header/templates_on.h"
#include "reim_source.c"
#include "../header/templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../header/templates_on.h"
#include "reim_source.c"
#include "../header/templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT
