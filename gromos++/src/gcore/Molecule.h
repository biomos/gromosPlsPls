// gcore_Molecule.h

#ifndef INCLUDED_GCORE_MOLECULE
#define INCLUDED_GCORE_MOLECULE

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

namespace gmath{
class Vec;
}

using gmath::Vec;

namespace gcore{

  class MoleculeTopology;

  class Molecule{
    MoleculeTopology *d_mt;
    std::vector<Vec*> d_pos;
    
    // not implemented
    Molecule();
    Molecule &operator=(const Molecule &);

  public:
    Molecule(const MoleculeTopology &mt);
    
    Molecule(const Molecule &);
    ~Molecule();
    
    // Methods
    Vec &pos(int i);
    
    // Accessors
    int numAtoms()const;
    const Vec &pos(int i)const;
    const MoleculeTopology &topology() const; 
    
  }; /* class Molecule */

  inline Vec &Molecule::pos(int i){
    assert(i < this->numAtoms());
    return *d_pos[i];
  }

  inline const Vec &Molecule::pos(int i)const{
    assert (i < this->numAtoms());
    return *d_pos[i];
  }
  
} /* Namespace */ 
#endif

