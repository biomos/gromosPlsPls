// gcore_Exclusion.h

#ifndef INCLUDED_GCORE_EXCLUSION
#define INCLUDED_GCORE_EXCLUSION

namespace gcore{

class Exclusion_i;
/**
 * Class Exclusion
 * Purpose: stores a set of atom numbers for e.g. exclusions / 1-4 neighbours
 *
 * Description:
 * The Exclusion class stores a set of atom numbers for e.g. exclusions or
 * 1-4 neighbours of an atom.
 *
 * @class Exclusion
 * @author R. Buergi
 * @ingroup gcore
 * @sa gcore::AtomTopology
 */
class Exclusion{
  Exclusion_i *d_this;
 public:
  /**
   * Exclusion constructor
   */
  Exclusion();
  /** 
   * Exclusion copy constructor
   * @param & Exclusion to be copied
   */
  Exclusion(const Exclusion&);
  /**
   * Exclusion deconstructor
   */
  ~Exclusion();

  // Methods
  /**
   * Member operator = copies one set of exclusions to the other
   */
  Exclusion &operator=(const Exclusion &);
  /**
   * Method to remove atom number i from the Exclusion
   * @param i Atom number to be removed
   */
  void erase(int i);
  // remove atom i from exclusion list

  /**
   * Method to add an atom to the Exclusion
   * @param i Atom number to be added
   */
  void insert(int i);
  // add Atom i to exclusion list

  /**
   * Accessor, returns the number of Exclusions
   */
  int size() const;
  // number of exclusions

  /**
   * Accessor, returns the atom number of the i-th exclusion as a const
   */
  int atom(int i) const;
  // get atom number being number i in exclusion list
  
  /**
   * check whether it contains an atom
   */
  bool contains(int i) const;
};

}

#endif

