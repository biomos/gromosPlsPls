// CGType.h
#ifndef INCLUDED_CGTYPE
#define INCLUDED_CGTYPE

namespace gcore{
  /**
   * Class CGType
   * Purpose: contains the coarse grain Lennard Jones parameters for an AtomPair
   *
   * Description:
   * The different Lennard Jones parameters (C12, C6) are 
   * contained in an CGType for every possible AtomPair (defined by their
   * Integer Atom Codes).
   *
   * @class CGType
   * @author N. Schmid
   * @ingroup gcore
   * @sa gcore::GromosForceField
   */
class CGType
{
  double d_c12, d_c6;
  public:
  /**
   * CGType constructor
   * @param c12 C12, C12 for normal interacting particles
   * @param c6 C6, C6 for normal interacting particles
   */
  CGType(double c12=0, double c6=0):
    d_c12(c12), d_c6(c6){}
  /**
   * CGType copy constructor
   * @param l CGType to be copied
   */
  CGType(const CGType &l): d_c12(l.d_c12), d_c6(l.d_c6) {}
  /**
   * Member operator = copies one CGType into the other
   */
  CGType &operator=(const CGType &l) {
    d_c12 = l.c12(); d_c6 = l.c6();
    return *this;
  }
  /**
   * CGType deconstructor
   */
  ~CGType(){}
  /**
   * Accessor, returns C12 (for normal interacting particles)
   */
  double c12()const{return d_c12;}
  /**
   * Accessor, returns C6 (for normal interacting particles)
   */
  double c6()const {return d_c6;}
};

}
#endif
