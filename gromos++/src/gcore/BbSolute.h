// gcore_BbSolute.h

#ifndef INCLUDED_GCORE_BBSOLUTE
#define INCLUDED_GCORE_BBSOLUTE

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{

  class BbSolute_i;
  class GromosForceField;
  class AtomTopology;
  class Exclusion;
  class Bond;
  class Angle;
  class Dihedral;
  class Improper;
  class BbBondIt;
  class BbAngleIt;
  class BbImpIt;
  class BbDihIt;
  /**
   * Class BbSolute
   * Purpose: defines a Building block for solute molecules (MTBUILDBLSOLUTE)
   *
   * Description:
   * Defines a Building block that is (part of) a solute molecule. Note 
   * that for the atoms there is a direct accessor, but that the bonds,
   * angles etc. are accessed via the iterators BbBondIt, BbAngleIt, 
   * BbDihIt and BbImpIt.
   *
   * @class BbSolute
   * @author C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::BbBondIt
   * @sa gcore::BbAngleIt
   * @sa gcore::BbDihIt
   * @sa gcore::BbImpIt
   * @sa gcore::BuildingBlock
   */
  class BbSolute{
    BbSolute_i *d_this;
    // This class contains all topological information
    /**
     * Bond Iterator for the BbSolute
     *
     * The BbSolute Bond iterator is used to loop over the Bonds in a
     * solute building block. 
     * It is constructed with the BbSolute as an argument. Use the ++ 
     * operator to move to the next Bond. The () operator returns the 
     * current Bond. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the bond list.
     * @author C. Oostenbrink
     */
    friend class BbBondIt;
    /**
     * Angle Iterator for the BbSolute
     *
     * The BbSolute Angle iterator is used to loop over the Angles in a
     * solute building block. 
     * It is constructed with the BbSolute as an argument. Use the ++ 
     * operator to move to the next Angle. The () operator returns the 
     * current Angle. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the Angle list.
     * @author C. Oostenbrink
     */
    friend class BbAngleIt;
    /**
     * Improper Iterator for the BbSolute
     *
     * The BbSolute Improper iterator is used to loop over the Impropers in a
     * solute building block. 
     * It is constructed with the BbSolute as an argument. Use the ++ 
     * operator to move to the next Improper. The () operator returns the 
     * current Improper. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the Improper list.
     * @author C. Oostenbrink
     */
    friend class BbImpIt;
    /**
     * Dihedral Iterator for the BbSolute
     *
     * The BbSolute Dihedral iterator is used to loop over the Dihedrals in a
     * solute building block. 
     * It is constructed with the BbSolute as an argument. Use the ++ 
     * operator to move to the next Dihedral. The () operator returns the 
     * current Dihedral. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the Dihedral list.
     * @author C. Oostenbrink
     */
    friend class BbDihIt;
    


  public:
    /**
     * BbSolute Constructor
     */
    BbSolute();
    /**
     * BbSolute copy constructor
     * @param & BbSolute to be copied
     */
    BbSolute(const BbSolute &);
    /**
     * BbSolute deconstructor
     */
    ~BbSolute();
    /**
     * Member operator = to copy one BbSolute into the other
     */
    BbSolute &operator=(const BbSolute &);
   
    // Methods
    /**
     * Member function to add an atom to the building block
     * @param a AtomTopology to be added; should be complete already
     */
    void addAtom(const AtomTopology &a);
    /**
     * Member function to add previous exclusions to the building block.
     * Previous exclusions are the exclusions that should be applied to the
     * previous building block
     * @param a Exclusions of type gcore::Exclusion
     */
    void addPexcl(const Exclusion &a);
    /**
     * Member function to add bonds to the building block
     * &param b Bond to be added; should be complete already
     */
    void addBond(const Bond &b);
    /**
     * Member function to add angles to the building block
     * &param b Angle to be added; should be complete already
     */
    void addAngle(const Angle &b);
    /**
     * Member function to add dihedral angles to the building block
     * &param b Dihedral to be added; should be complete already
     */
    void addDihedral(const Dihedral &b);
    /**
     * Member function to add improper dihedrals to the building block
     * @param b Improper to be added; should be complete already
     */
    void addImproper(const Improper &b);
    /**
     * Member function to set the name of the building block
     * @param s String containing the name
     */
    void setResName(const std::string &s);
    
    // Accessors
    /**
     * Accessor, returns the number of atoms in the building block
     */
    int numAtoms()const;
    /** 
     * Accessor, returns the number of preceding exclusions in the building
     * block
     */
    int numPexcl()const;
    /**
     * Accessor, returns an AtomTopology of atom i in the building block
     */
    const AtomTopology& atom(int i) const; 
    /** 
     * Accessor, returns the i-th set of preceding exclusions
     */
    const Exclusion& pexcl(int i) const;
    /**
     * Accessor, returns the name of the building block
     */
    const std::string &resName()const;
    
  }; /* class BbSolute */

  
  class BbBondIt_i;
  class BbAngleIt_i;
  class BbImpIt_i;
  class BbDihIt_i;

  class BbBondIt{
    BbBondIt_i *d_this;
    // not implemented
    BbBondIt();
    BbBondIt(const BbBondIt&);
    BbBondIt &operator=(const BbBondIt &);
  public:
    /**
     * BbBondIt constructor
     * @param mt is a solute building block
     */
    BbBondIt(const BbSolute &mt);
    /**
     * BeBondIt deconstructor
     */
    ~BbBondIt();
    /**
     * Member operator ++, moves to the next bond
     */
    void operator++();
    /**
     * Member operator (), returns the bond the iterator is currently 
     * pointing at
     */
    const Bond &operator()()const;
    /**
     * Member bool(), returns 1 as long as the iterator is not at the 
     * end of list of bonds in the building block
     */
    operator bool()const;
  };

  class BbAngleIt{
    BbAngleIt_i *d_this;
    // not implemented
    BbAngleIt();
    BbAngleIt(const BbAngleIt&);
    BbAngleIt &operator=(const BbAngleIt &);
  public:
    /**
     * BeAngleIt constructor
     * @param mt is an end-group building block
     */
    BbAngleIt(const BbSolute &mt);
    /**
     * BeAngleIt deconstructor
     */
    ~BbAngleIt();
    /**
     * Member operator ++, moves to the next angle
     */
    void operator++();
     /**
     * Member operator (), returns the angle the iterator is currently 
     * pointing at
     */
    const Angle &operator()()const;

     /**
     * Member bool(), returns 1 as long as the iterator is not at the 
     * end of list of angles in the building block
     */
   operator bool()const;
  };

  class BbImpIt{
    BbImpIt_i *d_this;
    // not implemented
    BbImpIt();
    BbImpIt(const BbImpIt&);
    BbImpIt &operator=(const BbImpIt &);
  public:
    /**
     * BbImpIt constructor
     * @param mt is an end-group building block
     */
    BbImpIt(const BbSolute &mt);
    /**
     * BbImpIt deconstructor
     */
    ~BbImpIt();
    /**
     * Member operator ++, moves to the next improper
     */
    void operator++();
    /**
     * Member operator (), returns the improper the iterator is currently 
     * pointing at
     */
    const Improper &operator()()const;
    /**
     * Member bool(), returns 1 as long as the iterator is not at the 
     * end of list of impropers in the building block
     */
    operator bool()const;
  };

  class BbDihIt{
    BbDihIt_i *d_this;
    // not implemented
    BbDihIt();
    BbDihIt(const BbDihIt&);
    BbDihIt &operator=(const BbDihIt &);
  public:
    /**
     * BbDihIt constructor
     * @param mt is an end-group building block
     */
    BbDihIt(const BbSolute &mt);
    /**
     * BbDihIt deconstructor
     */
    ~BbDihIt();
   /**
     * Member operator ++, moves to the next dihedral
     */ 
    void operator++();
    /**
     * Member operator (), returns the dihedral the iterator is currently 
     * pointing at
     */
    const Dihedral &operator()()const;
    /**
     * Member bool(), returns 1 as long as the iterator is not at the 
     * end of list of dihedrals in the building block
     */
    operator bool()const;
  };
  
} /* Namespace */ 
#endif



