// VirtualAtomType.h
#ifndef INCLUDED_VIRTUALATOMTYPE
#define INCLUDED_VIRTUALATOMTYPE

namespace gcore{

  /**
   * Class VirtualAtomType
   * Purpose: contains a virtual atom type
   *
   * Description:
   * Contains two distances that may be relevant to calculate the position of a virtual atom
   *
   * @class VirtualAtomType
   * @author R. Buergi, C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::GromosForceField
   */

class VirtualAtomType
{
  int d_code;
  double d_dis1;
  double d_dis2;
 public:
  /**
   * VirtualAtomType constructor
   * @param c    bond code
   * @param dis1 distance 1 (DISH)
   * @param dis2 distance 2 (DISC)
   */
  VirtualAtomType(int c=0, double dis1=0, double dis2=0): d_code(c), d_dis1(dis1), d_dis2(dis2){}
  /**
   * VirtualAtomType copyconstructor
   * @param va VirtualAtomType to be copied
   */
  VirtualAtomType(const VirtualAtomType& va): d_code(va.d_code), d_dis1(va.d_dis1), d_dis2(va.d_dis2){}
  /** 
   * Member operator=, assign distances of one virtual atom type to another
   */
  VirtualAtomType &operator=(const VirtualAtomType &va);
  /**
   * BondType deconstuctor
   */
  ~VirtualAtomType(){}
  /**
   * Accessor, returns the integer code
   */
  int code()const{return d_code;}
  /**
   * Accessor, returns the first distance
   */
  double dis1()const{return d_dis1;}
  /**
   * Accessor, returns the second distance
   */
  double dis2()const{return d_dis2;}
};

}
#endif



