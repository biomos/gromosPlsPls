// AngleType.h
#ifndef INCLUDED_ANGLETYPE
#define INCLUDED_ANGLETYPE

namespace gcore{

  /**
   * Class AngleType
   * Purpose: contains a gromos96 angle type
   *
   * Description:
   * Contains the optimum angle and force constant for a gromos96 bond
   * angle. The potential energy for an angle is defined as
   * @f[ V^{angle}=\frac{1}{2}K_{\theta_n}\left[\cos{\theta_n} - \cos{\theta_{0_n}}\right]^2@f]
   *
   * @class AngleType
   * @author R. Buergi
   * @sa gcore::Angle
   * @sa gcore::GromosForceField
   */
class AngleType
{
  double d_t0;
  double d_fc;
 public:
  /**
   * AngleType constructor
   * @param fc force constant (@f$K_{\theta_n}@f$)
   * @param l  optimum angle (@f$\theta_{0_n}@f$)
   */
  AngleType(double fc=0, double l=0): d_t0(l), d_fc(fc){}
  /**
   * AngleType copyconstructor
   * @param b AngleType to be copied
   */
  AngleType(const AngleType& b):d_t0(b.d_t0), d_fc(b.d_fc){}
  /** 
   * Member operator=, assign force constant and optimum angle of one
   * AngleType to the other
   */
  AngleType &operator=(const AngleType &b);
  /**
   * AngleType deconstuctor
   */
  ~AngleType(){}
  /**
   * Accessor, returns the optimum angle (@f$\theta_{0_n}@f$)
   */
  double t0()const{return d_t0;}
  /**
   * Accessor, returns the force constant (@f$K_{\theta_n}@f$)
   */
  double fc()const{return d_fc;}
};

}
#endif
