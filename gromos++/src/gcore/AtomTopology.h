// gcore_AtomTopology.h
#ifndef INCLUDED_GCORE_ATOMTOPOLOGY
#define INCLUDED_GCORE_ATOMTOPOLOGY

#ifndef INDLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{

class AtomTopology_i; 
class Exclusion; 

/**
 * Class AtomTopology
 * Purpose: contains all atomic properties for an atom
 *
 * Description:
 * All topological information that can be assigned to a single atom
 * is contained in an AtomTopology
 *
 * @class AtomTopology
 * @author R. Buergi
 * @sa gcore::Exclusion
 * @sa gcore::MoleculeTopology
 */
class AtomTopology {
  AtomTopology_i *d_this;

 public:
  // Constructors
  /**
   * AtomTopology constructor
   */
  AtomTopology();
  /**
   * AtomTopology copyconstructor
   * @param & AtomTopology to be copied
   */
  AtomTopology(const AtomTopology&);

  /**
   * AtomTopology deconstructor
   */
  ~AtomTopology();

  // Methods
  /**
   * Member operator = copies one AtomTopology into the other
   * @param & AtomTopology to be copied
   */
  AtomTopology &operator=(const AtomTopology&);
  /**
   * Member function to set the Integer Atom Code of the atom
   */
  void setIac(int);
  /**
   * Member function to set the ChargeGroup of the atom (0 or 1)
   */
  void setChargeGroup(int);
  /**
   * Member function to set the charge of the atom
   */
  void setCharge(double);
  /**
   * Member function to set the mass of the atom
   */
  void setMass(double);
  /**
   * Member function to set the name of the atom
   */
  void setName(const std::string&);
  /**
   * Member function to set the exclusions of the atom
   * @param e A set of exclusions of type gcore::Exclusion
   */
  void setExclusion(const Exclusion &e);
  /**
   * Member function to set the 1,4 neighbours of the atom
   * @param e A set of 1,4 neighbours of type gcore::Exclusion
   */
  void setExclusion14(const Exclusion &e);
  /**
   * Member function to set the Vanderwaals radius of the atom. This 
   * is not standard gromos96 topological information, but is used by 
   * some programs
   */
  void setradius(double);

  // Accessors
  /**
   * accessor, returns the Integer Atom Code of the atom
   */
  int iac()const;
  /**
   * accessor, returns the chargeGroup code of the atom (0 or 1)
   */
  int chargeGroup()const;
  /**
   * accessor, returns the charge of the atom
   */
  double charge()const;
  /**
   * accessor, returns the mass of the atom
   */
  double mass()const;
  /** 
   * accessor, returns the name of the atom
   */
  const std::string &name()const;
  /** 
   * accessor, returns the exclusions of the atom
   */
  const Exclusion &exclusion()const;
  /** 
   * accessor, returns the 1,4 neighbours of the atom
   */
  const Exclusion &exclusion14()const;
  /**
   * accessor, returns the Vanderwaals radius of the atom
   */
  double radius()const;
};

}
#endif

