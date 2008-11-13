// DihedralType.h
#ifndef INCLUDED_DIHEDRALTYPE
#define INCLUDED_DIHEDRALTYPE

namespace gcore {

/**
   * Class DihedralType<br>
   * Purpose: contains a gromos96 Dihedral angle type
   *
   * Description (still to be updated!):
   * Contains the phase (@f$\cos(\delta_n)@f$), multiplicity (@f$m_n@f$)
   * and force constant (@f$K_{\phi_n}@f$) for a gromos96 Dihedral angle.
   * The potential energy for a trigonomic dihedral angle is defined as
   * @f[ V^{trig}=K_{\phi_n}\left[1+\cos(\delta_n)\cos(m_n\phi_n)\right]
   * @f]
   *
   * @class DihedralType
   * @author R. Buergi, D. Geerke
   * @ingroup gcore
   * @sa gcore::Dihedral
   * @sa gcore::GromosForceField
   */
  class DihedralType {
    int d_code;
    double d_fc;
    double d_pd;
    double d_pdl;
    int d_np;
  public:

    /**
     * DihedralType constructor
     * @param c  integer code of the dihedral type
     * @param fc Force constant     (@f$K_{\phi_n}@f$)
     * @param pd phase              (@f$\cos(\delta_n)@f$)
     * @param pdl phase-shift angle
     * @param np multiplicity       (@f$m_n@f$)
     */
    DihedralType(int c = 0, double fc = 0, double pd = 0, double pdl = 0, int np = 0) :
    d_code(c), d_fc(fc), d_pd(pd), d_pdl(pdl), d_np(np) {
    }

    /**
     * DihedralType copy constructor
     * @param b Dihedral to be copied
     */
    DihedralType(const DihedralType& b) : d_code(b.d_code), d_fc(b.d_fc),
            d_pd(b.d_pd), d_pdl(b.d_pdl), d_np(b.d_np) {
    }
    /**
     * Member operator = copies two DihedralTypes
     */
    DihedralType & operator=(const DihedralType &b);

    /** 
     * Accessor, returns the inetger code of the dihedral type
     */
    int code()const {
      return d_code;
    }
    /**
     * Accessor, returns the phase (@f$\cos(\delta_n)@f$)
     */
    double pd()const {
      return d_pd;
    }

    /**
     * Accessor, returns the multiplicity (@f$m_n@f$)
     */
    int np()const {
      return d_np;
    }

    /**
     * Accessor, returns the Force constant (@f$K_{\phi_n}@f$)
     */
    double fc()const {
      return d_fc;
    }

    /**
     * Accessor, returns the phase-shift angle in degrees
     */
    double pdl()const {
      return d_pdl;
    }
  };

}
#endif
