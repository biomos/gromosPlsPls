// BondType.h
#ifndef INCLUDED_BONDTYPE
#define INCLUDED_BONDTYPE

namespace gcore{

  /**
   * Class BondType
   * Purpose: contains a gromos96 bond type
   *
   * Description:
   * Contains the optimum bond length (@f$b_{0_n}@f$) and force constant
   * (@f$K_{b_n}@f$) for a gromos96 bond. The potential energy for a bond 
   * is defined as
   * @f[ V^{bond}=\frac{1}{4}K_{b_n}\left[b_n^2 - b_{0_n}^2\right]^2@f]
   *
   * @class BondType
   * @author R. Buergi
   * @sa gcore::Bond
   * @sa gcore::GromosForceField
   */

class BondType
{
  double d_fc;
  double d_b0;
 public:
  /**
   * BondType constructor
   * @param fc force constant (@f$K_{b_n}@f$)
   * @param l  optimum angle (@f$b_{0_n}@f$)
   */
  BondType(double fc=0, double l=0): d_fc(fc), d_b0(l){}
  /**
   * BondType copyconstructor
   * @param b BondType to be copied
   */
  BondType(const BondType& b):d_fc(b.d_fc), d_b0(b.d_b0){}
  /** 
   * Member operator=, assign force constant and optimum bond lenght of one
   * BondType to the other
   */
  BondType &operator=(const BondType &b);
  /**
   * BondType deconstuctor
   */
  ~BondType(){}
  /**
   * Accessor, returns the optimum bond length (@f$b_{0_n}@f$)
   */
  double b0()const{return d_b0;}
  /**
   * Accessor, returns the force constant (@f$K_{b_n}@f$)
   */
  double fc()const{return d_fc;}
};

}
#endif



