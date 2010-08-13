// fit_AtomDistances.cc
//includes explicit calls to gsl now

#include <cassert>

#include <iostream>
#include <vector>
#include <sstream>
#include <set>

#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Matrix.h"
#include "../gmath/Vec.h"
#include "PositionUtils.h"
#include "../utils/AtomSpecifier.h"

#include "AtomDistances.h"

using gmath::Matrix;
using gmath::Vec;
using namespace fit;
using namespace std;



double AtomDistances::dist(vector<Vec> const &ref, 
			       vector<Vec> const &sys)const
{
  double dist2=0.0;
  for(size_t i=0; i< ref.size(); ++i){
    for(size_t j=i+1; j < ref.size(); ++j){
      for(int b=0;b<3;++b){
	dist2 += 
	  ((ref[i][b] - ref[j][b]) - (sys[i][b] - sys[j][b])) *
	  ((ref[i][b] - ref[j][b]) - (sys[i][b] - sys[j][b]));
      }
    }
  }
  int N = ref.size() * (ref.size() - 1) / 2;
  
  return sqrt(dist2/N);
}
