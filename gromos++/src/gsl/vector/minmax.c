#include "../../../config.h"
#include <stdlib.h>
#include "gsl_vector.h"

#define BASE_LONG_DOUBLE
#include "../header/templates_on.h"
#include "minmax_source.c"
#include "../header/templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../header/templates_on.h"
#include "minmax_source.c"
#include "../header/templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../header/templates_on.h"
#include "minmax_source.c"
#include "../header/templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../header/templates_on.h"
#include "minmax_source.c"
#include "../header/templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../header/templates_on.h"
#include "minmax_source.c"
#include "../header/templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../header/templates_on.h"
#include "minmax_source.c"
#include "../header/templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../header/templates_on.h"
#include "minmax_source.c"
#include "../header/templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../header/templates_on.h"
#include "minmax_source.c"
#include "../header/templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../header/templates_on.h"
#include "minmax_source.c"
#include "../header/templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../header/templates_on.h"
#include "minmax_source.c"
#include "../header/templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../header/templates_on.h"
#include "minmax_source.c"
#include "../header/templates_off.h"
#undef  BASE_CHAR


