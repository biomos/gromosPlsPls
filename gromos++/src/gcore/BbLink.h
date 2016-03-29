// gcore_BbSolute.h
#ifndef INCLUDED_BBLINK
#define INCLUDED_BBLINK

namespace gcore{

  class GromosForceField;
  class AtomTopology;
  class Bond;
  class Angle;
  class Dihedral;
  class Improper;
  class BbSolute;
  /**
   * Class BbLINK
   * Purpose: defines a Building block for linking 
   *
   *
   * @class BbLink
   * @author C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::MoleculeTopology
   * @sa gcore::BuildingBlock
   */
  class BbLink: public gcore::BbSolute{
    std::vector<int> d_linkres;
  public:
    /**
     * BbSolute Constructor
     */
    BbLink(){setRep(0);};
    /**
     * BbSolute copy constructor
     */
    BbLink(const BbLink &);
    /**
     * BbSolute deconstructor
     */
    ~BbLink(){};
    /**
     * Member operator = to copy one BbSolute into the other
     */
    BbLink &operator=(const BbLink &);

    // Methods
    /**
     * Member function to store the residue identifier for a given atom
     * @param a atom
     * @param i residue
     */
    void setLinkRes(const int a, const int i);

    // Accessors
    /** 
     * Accessor, returns the residue identifier for a given atom a
     */
    int linkRes(int a)const;

  }; /* class BbLink */
} /* Namespace */ 



#endif
