// utils_Neighbours.h

// Class that contains bonded neighbours of atom i

#ifndef INCLUDED_UTILS_NEIGHBOURS
#define INCLUDED_UTILS_NEIGHBOURS

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

namespace gcore{
  class Simulation;
  class System;
  class Molecule;
  class MoleculeTopology;
}

namespace utils{
class Neighbours: public std::vector<int>{
  // not implemented
  Neighbours (const Neighbours&);
  Neighbours();
  Neighbours &operator=(const Neighbours&);
 public:
  Neighbours(const gcore::System &sys, int mol, int i);
  Neighbours(const gcore::Molecule &mol, int k);
  Neighbours(const gcore::MoleculeTopology &mol, int k);
};
}
#endif
