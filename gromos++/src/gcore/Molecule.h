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

/**
 * Class Molecule
 * Purpose: Contains coordinates and topology of a Molecule
 *
 * Description:
 * The gromos++ Molecule is a specific entity. It is defined as a set of 
 * atoms that are connected by bonds. The class Molecule contains both the
 * coordinates of the molecule and the topological information
 *
 * @class Molecule
 * @author R. Buergi
 * @sa gcore::MoleculeTopology
 * @sa gcore::System
 */
namespace gcore{

  class MoleculeTopology;

  class Molecule{
    MoleculeTopology *d_mt;
    std::vector<Vec*> d_pos;
    
    // not implemented
    Molecule();
    Molecule &operator=(const Molecule &);

  public:
    /**
     * Molecule constructor
     * @param mt is a MoleculeTopology, a new molecule (with no coordinates)
     *           will be constructed based on this MoleculeTopology
     */
    Molecule(const MoleculeTopology &mt);
    /**
     * Molecule copy constructor
     * @param & the molecule to be copied
     */
    Molecule(const Molecule &);
    /**
     * Molecule deconstructor
     */
    ~Molecule();
    
    // Methods
    /**
     * Accessor, returns the a pointer to a vector with the coordinates 
     * of the i-th atom in the Molecule
     * @return A gmath::Vec of three coordinates
     */
    Vec &pos(int i);
    
    // Accessors
    /**
     * Accessor, returns the number of atoms in the Molecule
     */
    int numAtoms()const;
    /**
     * Accessor, returns the coordinates of the i-th atom in the Molecule 
     * as a const
     * @return A gmath::Vec of three coordinates
     */
    const Vec &pos(int i)const;
    /**
     * Accessor, returns the MoleculeTopology of the Molecule
     */
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

