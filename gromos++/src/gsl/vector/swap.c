#include "../../../config.h"
#include "../err/gsl_errno.h"
#include "gsl_vector.h"

#define BASE_GSL_COMPLEX_LONG
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../header/../header/templates_on.h"
#include "swap_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_CHAR
