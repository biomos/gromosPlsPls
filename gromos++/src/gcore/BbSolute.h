// gcore_BbSolute.h
#ifndef INCLUDED_GCORE_BBSOLUTE
#define INCLUDED_GCORE_BBSOLUTE

#ifndef INCLUDED_GCORE_MOLECULETOPOLOGY
#include "MoleculeTopology.h"
#define INCLUDED_GCORE_MOLECULETOPOLOGY
#endif
#ifndef INCLUDED_GCORE_EXCLUSION
#include "Exclusion.h"
#define INCLUDED_GCORE_EXCLUSION
#endif
#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif
#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

namespace gcore{

  class GromosForceField;
  class MoleculeTopology;
  class AtomTopology;
  class Exclusion;
  class Bond;
  class Angle;
  class Dihedral;
  class Improper;
  /**
   * Class BbSolute
   * Purpose: defines a Building block for solute molecules (MTBUILDBLSOLUTE)
   *
   * Description:
   * Defines a Building block that is (part of) a solute molecule. Note 
   * that for the atoms there is a direct accessor, but that the bonds,
   * angles etc. are accessed via the iterators, just as in the 
   * MoleculeTopology.
   *
   * @class BbSolute
   * @author C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::MoleculeTopology
   * @sa gcore::BuildingBlock
   */
  class BbSolute: public gcore::MoleculeTopology{
    std::vector<Exclusion> d_pexcl;
    int d_rep;
    
  public:
    /**
     * BbSolute Constructor
     */
    BbSolute(){d_rep=0;};
    /**
     * BbSolute copy constructor
     */
    BbSolute(const BbSolute &);
    
    /**
     * BbSolute deconstructor
     */
    ~BbSolute(){};
    /**
     * Member operator = to copy one BbSolute into the other
     */
    BbSolute &operator=(const BbSolute &);

    // Methods
    /**
     * Member function to add previous exclusions to the building block.
     * Previous exclusions are the exclusions that should be applied to the
     * previous building block
     * @param a Exclusions of type gcore::Exclusion
     */
    void addPexcl(const Exclusion &a);
    /**
     * Member function to set the name of the building block
     * @param s String containing the name
     */
    void setResName(const std::string &s);
    /**
     * Member function to set the number of atoms that replace atoms in
     * linked building blocks. If i<0 the first i atoms replace the last
     * i atoms of the previous building block. If i>0 the last i atoms
     * replace the first i atoms of the following building block.<br>
     * see also the documentation on the <a
     * href=http://wwi.igc.ethz.ch/gromos++/maketop.html>program maketop</a>.
     */
    void setRep(const int i);

    // Accessors
    /** 
     * Accessor, returns the number of preceding exclusions in the building
     * block
     */
    int numPexcl()const;
    /** 
     * Accessor, returns the i-th set of preceding exclusions
     */
    const Exclusion& pexcl(int i) const;
    /**
     * Accessor, returns the name of the building block
     */
    const std::string &resName()const;
    /**
     * Accessor, returns the number of atoms to be replaced. If i<0 the
     * first i atoms replace the last i atoms of the previous building
     * block. If i>0 the last i atoms replace the first i atoms of the
     * following building block.<br>
     * see also the documentation on the <a
     * href=http://wwi.igc.ethz.ch/gromos++/maketop.html>program maketop</a>.
     */
    const int rep()const;

  }; /* class BbSolute */
} /* Namespace */ 
#endif



