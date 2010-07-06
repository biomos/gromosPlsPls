// LJType.h
#ifndef INCLUDED_LJTYPE
#define INCLUDED_LJTYPE

namespace gcore{
  /**
   * Class LJType
   * Purpose: contains the Lennard Jones parameters for an AtomPair
   *
   * Description:
   * The different Lennard Jones parameters (C12, C6, CS12 and CS6) are 
   * contained in an LJType for every possible AtomPair (defined by their
   * Integer Atom Codes).
   *
   * @class LJType
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::GromosForceField
   */
class LJType
{
  double d_c12, d_c6, d_cs12, d_cs6;
  public:
  /**
   * LJType constructor
   * @param d1 C12, C12 for normal interacting particles
   * @param d2 C6, C6 for normal interacting particles
   * @param d3 CS12, C12 for 1,4 neighbours
   * @param d4 CS6, C6 for 1,4 neighbours
   */
  LJType(double d1=0, double d2=0, double d3=0, double d4=0):
    d_c12(d1), d_c6(d2), d_cs12(d3), d_cs6(d4){}
  /**
   * LJType copy constructor
   * @param l LJType to be copied
   */
  LJType(const LJType &l): d_c12(l.d_c12), d_c6(l.d_c6), 
    d_cs12(l.d_cs12), d_cs6(l.d_cs6){}
  /**
   * Member operator = copies one LJType into the other
   */
  LJType &operator=(const LJType &l);
  /**
   * LJType deconstructor
   */
  ~LJType(){}
  /**
   * Accessor, returns C12 (for normal interacting particles)
   */
  double c12()const{return d_c12;}
  /**
   * Accessor, returns C6 (for normal interacting particles)
   */
  double c6()const {return d_c6;}
  /**
   * Accessor, returns CS12 (C12 for 1,4 neighbours)
   */
  double cs12()const {return d_cs12;}
  /**
   * Accessor, returns CS6 (C6 for 1,4 neighbours)
   */
  double cs6()const {return d_cs6;}
   /**
   * Accessor, returns C12 (for normal interacting particles)
   */
  double & c12(){return d_c12;}
  /**
   * Accessor, returns C6 (for normal interacting particles)
   */
  double & c6() {return d_c6;}
  /**
   * Accessor, returns CS12 (C12 for 1,4 neighbours)
   */
  double & cs12() {return d_cs12;}
  /**
   * Accessor, returns CS6 (C6 for 1,4 neighbours)
   */
  double & cs6() {return d_cs6;}
};

}
#endif
