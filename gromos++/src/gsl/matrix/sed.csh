#!/bin/csh

foreach i (*.h)
 sed 's/<gsl\/gsl_errno.h>/"..\/err\/gsl_errno.h"/g' $i >! tmp
 mv tmp $i
 sed 's/<gsl\/gsl_complex.h>/"..\/complex\/gsl_complex.h"/g' $i >! tmp
 mv tmp $i
 sed 's/<gsl\/gsl_block.h>/"..\/block\/gsl_block.h"/g' $i >! tmp
 mv tmp $i
 sed 's/templates_on.h/..\/header\/templates_on.h/g' $i >! tmp
 mv tmp $i
 sed 's/templates_off.h/..\/header\/templates_off.h/g' $i >! tmp
 mv tmp $i
 sed 's/<config.h>/"..\/..\/..\/config.h"/g' $i >! tmp
 mv tmp $i
 sed 's/<gsl\/gsl_matrix.h>/"..\/matrix\/gsl_matrix.h"/g' $i >! tmp
 mv tmp $i
 sed 's/<gsl\/gsl_vector.h>/"..\/vector\/gsl_vector.h"/g' $i >! tmp
 mv tmp $i
 sed 's/<gsl\/gsl_math.h>/"..\/header\/gsl_math.h"/g' $i >! tmp
 mv tmp $i
 sed 's/<gsl\/gsl_permutation.h>/"..\/permutation\/gsl_permutation.h"/g' $i >! tmp
 mv tmp $i
end

foreach i (*.c)
 sed 's/<gsl\/gsl_errno.h>/"..\/err\/gsl_errno.h"/g' $i >! tmp
 mv tmp $i
 sed 's/<gsl\/gsl_complex.h>/"..\/complex\/gsl_complex.h"/g' $i >! tmp
 mv tmp $i
 sed 's/<gsl\/gsl_block.h>/"..\/block\/gsl_block.h"/g' $i >! tmp
 mv tmp $i
  sed 's/templates_on.h/..\/header\/templates_on.h/g' $i >! tmp
 mv tmp $i
 sed 's/templates_off.h/..\/header\/templates_off.h/g' $i >! tmp
 mv tmp $i
  sed 's/<config.h>/"..\/..\/..\/config.h"/g' $i >! tmp
 mv tmp $i
  sed 's/<gsl\/gsl_matrix.h>/"..\/matrix\/gsl_matrix.h"/g' $i >! tmp
 mv tmp $i
 sed 's/<gsl\/gsl_vector.h>/"..\/vector\/gsl_vector.h"/g' $i >! tmp
 mv tmp $i
 sed 's/<gsl\/gsl_math.h>/"..\/header\/gsl_math.h"/g' $i >! tmp
 mv tmp $i
 sed 's/<gsl\/gsl_permutation.h>/"..\/permutation\/gsl_permutation.h"/g' $i >! tmp
 mv tmp $i
end

