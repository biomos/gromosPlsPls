// gcore_BbEnd.h

#ifndef INCLUDED_GCORE_BBEND
#define INCLUDED_GCORE_BBEND

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{

  class BbEnd_i;
  class AtomTopology;
  class Exclusion;
  class Bond;
  class Angle;
  class Dihedral;
  class Improper;
  class BeBondIt;
  class BeAngleIt;
  class BeImpIt;
  class BeDihIt;

  /**
   * Class BbEnd
   * Purpose: defines a Building block that is an end-group of a chain
   *
   * Description:
   * Defines a Building block that is an end-group of a chain. Very similar
   * the class BbSolute. Not that for the atoms there is a direct accessor,
   * but the bonds etc. are accessed via the iterators, BeBondIt, BeAngleIt,
   * BeDihIt and BeImpIt
   *
   * @class BbEnd
   * @author C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::BeBondIt
   * @sa gcore::BeAngleIt
   * @sa gcore::BeDihIt
   * @sa gcore::BeImpIt
   * @sa gcore::BuildingBlock
   */

  class BbEnd{
    BbEnd_i *d_this;
    // This class contains all topological information

    /**
     * Bond Iterator for the BbEnd
     *
     * The BbEnd bond iterator is used to loop over the bonds in an 
     * end-group building block. 
     * It is constructed with the BbEnd as an argument. Use the ++ operator 
     * to move to the next Bond. The () operator returns the current Bond. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the Bond list.
     * @author C. Oostenbrink
     */
    friend class BeBondIt;
    /**
     * Angle Iterator for the BbEnd
     *
     * The BbEnd angle iterator is used to loop over the angles in an 
     * end-group building block. 
     * It is constructed with the BbEnd as an argument. Use the ++ operator 
     * to move to the next Angle. The () operator returns the current Angle. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the Angle list.
     * @author C. Oostenbrink
     */
    friend class BeAngleIt;
    /**
     * Improper Iterator for the BbEnd
     *
     * The BbEnd Improper iterator is used to loop over the Impropers in an 
     * end-group building block. 
     * It is constructed with the BbEnd as an argument. Use the ++ operator 
     * to move to the next Improper. The () operator returns the current 
     * Improper. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the Improper list.
     * @author C. Oostenbrink
     */
    friend class BeImpIt;
    /**
     * Dihedral Iterator for the BbEnd
     *
     * The BbEnd Dihedral iterator is used to loop over the Dihedrals in an 
     * end-group building block. 
     * It is constructed with the BbEnd as an argument. Use the ++ operator 
     * to move to the next Dihedral. The () operator returns the current 
     * Dihedral. 
     * This can also be used as a boolean: the bool() returns 1 as long as 
     * the iterator is not at the end of the Dihedral list.
     * @author C. Oostenbrink
     */
    friend class BeDihIt;
    
  public:
    /**
     * BbEnd constructor
     */
    BbEnd();
    /**
     * BbEnd copy constructor
     * @param & BbEnd to be copied
     */
    BbEnd(const BbEnd &);
    /**
     * BbEnd deconstructor
     */
    ~BbEnd();
    /**
     * Member operator =, copies one BbEnd into the other
     */
    BbEnd &operator=(const BbEnd &);
    
    // Methods
    /**
     * Member function to add an atom to the building block
     * @param a AtomTopology to be added; should already be complete
     */
    void addAtom(const AtomTopology &a);
    /**
     * Member function to add a bond to the building block
     * @param b Bond to be added; should be complete
     */
    void addBond(const Bond &b);
    /** 
     * Member function to add an angle to the building block
     * @param b Angle to be added; should be complete
     */
    void addAngle(const Angle &b);
    /**
     * Member function to add a dihedral angle to the building block
     * @param b Dihedral to be added; should be complete
     */
    void addDihedral(const Dihedral &b);
    /**
     * Member function to add an improper dihedral to the building block
     * @param b Improper to be added; should be complete
     */
    void addImproper(const Improper &b);
    /**
     * Member function to set the name of the building block
     * @param s String with the name
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
     * Accessor, returns the number of atoms in the building block
     */
    int numAtoms()const;
    /**
     * Accessor, returns atom i of the building block
     * @param i atom number
     * @return An AtomTopology
     */
    const AtomTopology& atom(int i) const; 
    /**
     * Accessor returns the name of the building block
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
    
  }; /* class BbEnd */

  
  class BeBondIt_i;
  class BeAngleIt_i;
  class BeImpIt_i;
  class BeDihIt_i;

  class BeBondIt{
    BeBondIt_i *d_this;
    BeBondIt();
    BeBondIt(const BeBondIt&);
    BeBondIt &operator=(const BeBondIt &);
  public:
    BeBondIt(const BbEnd &mt);
    ~BeBondIt();
    void operator++();
    const Bond &operator()()const;
    operator bool()const;
  };

  class BeAngleIt{
    BeAngleIt_i *d_this;
    // not implemented
    BeAngleIt();
    BeAngleIt(const BeAngleIt&);
    BeAngleIt &operator=(const BeAngleIt &);
  public:
    BeAngleIt(const BbEnd &mt);
    ~BeAngleIt();
    void operator++();
    const Angle &operator()()const;
    operator bool()const;
  };

  class BeImpIt{
    BeImpIt_i *d_this;
    // not implemented
    BeImpIt();
    BeImpIt(const BeImpIt&);
    BeImpIt &operator=(const BeImpIt &);
  public:
    BeImpIt(const BbEnd &mt);
    ~BeImpIt();
    void operator++();
    const Improper &operator()()const;
    operator bool()const;
  };

  class BeDihIt{
    BeDihIt_i *d_this;
    // not implemented
    BeDihIt();
    BeDihIt(const BeDihIt&);
    BeDihIt &operator=(const BeDihIt &);
  public:
    BeDihIt(const BbEnd &mt);
    ~BeDihIt();
    void operator++();
    const Dihedral &operator()()const;
    operator bool()const;
  };
  
} /* Namespace */ 
#endif



