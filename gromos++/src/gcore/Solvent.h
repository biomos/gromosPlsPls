// gcore_Solvent.h

#ifndef INCLUDED_GCORE_SOLVENT
#define INCLUDED_GCORE_SOLVENT

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

namespace gmath{
class Vec;
}

using gmath::Vec;

namespace gcore{

  class SolventTopology;

  class Solvent{
    SolventTopology *d_mt;
    std::vector<Vec*> d_pos;
    int d_numCoords;
    
    // not implemented
    Solvent();
    //Solvent &operator=(const Solvent &);

  public:
    Solvent(const SolventTopology &mt);
    
    Solvent(const Solvent &);
    ~Solvent();
    
    // Methods
    Vec &pos(int i);
    void addCoord(Vec v);
    void setnumCoords(int i);
    Solvent &operator=(const Solvent &s);
    
    // Accessors
    int numCoords()const;
    
    const Vec &pos(int i)const;
    const SolventTopology &topology() const; 
    
  }; /* class Solvent */

  inline Vec &Solvent::pos(int i){
    assert(i < (this->numCoords()));
    return *d_pos[i];
  }

  inline const Vec &Solvent::pos(int i)const{
    assert (i < (this->numCoords()));
    return *d_pos[i];
  }
  
} /* Namespace */ 
#endif

