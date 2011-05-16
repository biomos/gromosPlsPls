// ImproperType.h
#ifndef INCLUDED_IMPROPERTYPE
#define INCLUDED_IMPROPERTYPE

namespace gcore{
  /**
   * Class ImproperType
   * Purpose: contains a gromos96 Improper dihedral angle type
   *
   * Description:
   * Contains the optimum angle value (@f$\xi_{0_n}@f$) and the force 
   * constant (@f$K_{\xi_n}@f$) for a gromos96 Improper angle.
   * The potential energy for a harmonic dihedral angle (Improper) is 
   * defined as
   * @f[ V^{har}=\frac{1}{2}K_{\xi_n}\left[\xi_n - \xi_{0_n}\right]^2
   * @f]
   *
   * @class ImproperType
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::Improper
   * @sa gcore::GromosForceField
   */
class ImproperType{
  int d_code;
  double d_q0;
  double d_fc;
 public:
  /**
   * ImproperType constructor
   * @param c  the integer code of the improper dihedral
   * @param fc The force constant (@f$K_{\xi_n}@f$)
   * @param l  The optimum angle value (@f$\xi_{0_n}@f$)
   */
  ImproperType(int c = 0,double fc=0, double l=0): d_code(c), d_q0(l),
          d_fc(fc){}
  /**
   * ImproperType copy constructor
   */
  ImproperType(const ImproperType& b): d_code(b.d_code), d_q0(b.d_q0),
          d_fc(b.d_fc){}
  /**
   * ImproperType deconstructor
   */
  ~ImproperType(){}
  /**
   * Member operator = copies one ImproperType into the other
   */
  ImproperType &operator=(const ImproperType &);
  /**
   * Accessor, returns the integer code
   */
  int code()const{return d_code;}
  /**
   * Accessor, returns the optimum angle value of the Improper 
   * (@f$\xi_{0_n}@f$)
   */
  double q0()const{return d_q0;}
  /** 
   * Accessor, returns the force constant of the Improper (@f$K_{\xi_n}@f$)
   */
  double fc()const{return d_fc;}
};

}
#endif
