#ifndef INCLUDED_UTILS_GCH
#define INCLUDED_UTILS_GCH
#endif
#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif
#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif
#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#define INCLUDED_GROMOS_EXCEPTION
#endif

#include <vector>
#include "../gcore/System.h"

namespace gcore{
  class System;
  class GromosForcefield;
  class Bond;
  class Angle;
}

namespace gmath{
  class Vec;
}

namespace utils
{
  /**
  * fill the hydrogen and non-hydrogen neighbours of atom a of molecule m
  * into vectors h and nh
  */
  void get_h_nh_neighbours(gcore::System &sys, gcore::GromosForceField &gff, int m, int a, std::vector<int> &h, std::vector<int> &nh);
   
  /**
  * determine the geometry of an atom based on the number of hydrogen
  * and non-hydrogen neighbours
  */ 
  int get_geometry(int numH, int numNH);

  /**
  * recalculate optimal positions for hydrogens based on the 
  * @param m molecule number
  * @param a atom number
  * @param h vector of hydrogen neighbours of a
  * @param nh vector of non-hydrogen neighbours of a
  * @param geom geometry type, see @ref gch 
  * @param eps tolerated deviation from optimal bond length
  */
  int generate_hcoordinates(gcore::System &sys, gcore::GromosForceField &gff, int m, int a, std::vector<int> &h, std::vector<int> &nh, int geom, double eps);  
  int generate_hcoordinates(gcore::System &sys, gcore::GromosForceField &gff, int m, int a, double eps);

}
