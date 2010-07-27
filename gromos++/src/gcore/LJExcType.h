// LJExcType.h
#ifndef INCLUDED_LJEXCTYPE
#define INCLUDED_LJEXCTYPE

namespace gcore{
  /**
   * Class LJExcType
   * Purpose: contains the Lennard Jones exception parameters of an AtomPair
   *
   * Description:
   * The different Lennard Jones exception parameters (C12, C6) are
   * contained in an LJExcType for the exceptions indicated in the LJEXCEPTIONS
   * block of the topology.
   *
   * @class LJExcType
   * @author A.P. Eichenberger
   * @ingroup gcore
   * @sa gcore::GromosForceField
   */
class LJExcType
{
  double d_c12, d_c6;
  public:
  /**
   * LJExcType constructor
   * @param d1 C12 parameter of the exceptioned pair
   * @param d2 C6 parameter of the exceptioned pair
   */
  LJExcType(double d1=0, double d2=0):
    d_c12(d1), d_c6(d2) {}
  /**
   * LJExcType copy constructor
   * @param l LJExcType to be copied
   */
  LJExcType(const LJExcType &l): d_c12(l.d_c12), d_c6(l.d_c6) {}
  /**
   * Member operator = copies one LJExcType into the other
   */
  LJExcType &operator=(const LJExcType &l);
  /**
   * LJExcType deconstructor
   */
  ~LJExcType(){}
  /**
   * Accessor, returns C12 (for normal interacting particles)
   */
  double c12()const{return d_c12;}
  /**
   * Accessor, returns C6 (for normal interacting particles)
   */
  double c6()const {return d_c6;}

   /**
   * Accessor, returns C12 (for normal interacting particles)
   */
  double & c12(){return d_c12;}
  /**
   * Accessor, returns C6 (for normal interacting particles)
   */
  double & c6() {return d_c6;}
};

}
#endif
