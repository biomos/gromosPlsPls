/* vector/vector.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "../../../config.h"
#include "../err/gsl_errno.h"
#include "gsl_vector.h"

/* turn off range checking at runtime if zero */
int gsl_check_range = 1;

#define BASE_GSL_COMPLEX_LONG
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../header/../header/templates_on.h"
#include "vector_source.c"
#include "../header/../header/templates_off.h"
#undef  BASE_CHAR
