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

  /** 
   * Class Solvent
   * Purpose: Contains coordinates and topological information of the solvent
   *
   * Description:
   * The gromos++ Solvent is a special kind of molecule. The class contains
   * topological information for the solvent (one) as well as the coordinates
   * for ALL solvent molecules in the system (many molecules). This of
   * course because it is not very useful to store 15000 copies of the same 
   * topology.
   *
   * @class Solvent
   * @author M. Kastenholz and C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::SolventTopology
   * @sa gcore::System
   */
  class Solvent{
    SolventTopology *d_mt;
    std::vector<Vec*> d_pos;
    int d_numCoords;
    
    // not implemented
    Solvent();
    //Solvent &operator=(const Solvent &);

  public:
    /**
     * Solvent Constructor
     * @param mt a SolventTopology, this creates a solvent without coordinates
     */
    Solvent(const SolventTopology &mt);
    /**
     * Solvent Copy Constructor
     * @param & Solvent to be copied
     */
    Solvent(const Solvent &);
    /**
     * Solvent deconstructor
     */
    ~Solvent();
    
    // Methods
    /**
     * Accessor, returns a pointer to a vector with the coordinates of 
     * solvent atom i. 
     * @param i The solvent atom of which you want the coordinates. i simply
     *          goes over all Solvent atoms (all atoms of all molecules)
     * @return A gmath::Vec containing three coordinates
     */
    Vec &pos(int i);
    /**
     * Method to add an atom to the Solvent.
     *
     * The user is responsible for adding complete molecules, that is adding
     * all atoms of solvent molecule to the Solvent (and in the correct order)
     * @param v A gmath::Vec containing three coordinates
     */
    void addCoord(Vec v);
    /**
     * Method to rescale the number of atoms in the Solvent
     *
     * This method allows you to set the total number of solvent atoms we have
     * If i < numCoords then the memory containing the coordinates for all 
     * atoms >i is released.
     * @param i The new total number of Solvent atoms in the class
     */
    void setnumCoords(int i);
    /**
     * Member operator = copies one Solvent into the other
     */
    Solvent &operator=(const Solvent &s);
    
    // Accessors
    /**
     * Accessor, returns the number of atoms for which coordinates are 
     * stored in the Solvent class. So this is the total number of solvent
     * atoms the class knows about (=number of atoms per solvent molecules
     * * number of solvent molecules)
     */
    int numCoords()const;
    /**
     * Accessor, returns a pointer to a vector with the coordinates of 
     * solvent atom i as a const
     * @param i The solvent atom of which you want the coordinates. i simply
     *          goes over all Solvent atoms (all atoms of all molecules)
     * @return A gmath::Vec containing three coordinates
     */
    const Vec &pos(int i)const;
    /**
     * Accessor, returns a SolventTopology containing the topological 
     * information for one (1) solvent molecule
     */
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

