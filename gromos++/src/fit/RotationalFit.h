// fit_RotationalFit.h

#ifndef INCLUDED_FIT_ROTATIONALFIT
#define INCLUDED_FIT_ROTATIONALFIT


#include "../gromos/Exception.h"

namespace gcore{
  class System;
}

namespace fit{
  class RotationalFit_i;
  class Reference;
  /**
   * Class RotationalFit
   * This class performs a rotational fit on the system
   *
   * A least squares fitting of one system is performed relative to a
   * reference system. The atoms that are taken into account are defined
   * by the Reference class
   *
   * @class RotationalFit
   * @author R. Buergi
   * @ingroup fit
   * @sa fit::TranslationalFit
   * @sa fit::PositionUtils
   * @sa utils::Rmsd
   */
  class RotationalFit{
    Reference *d_ref;

    // not implemented
    RotationalFit();
    RotationalFit(const RotationalFit &);
    RotationalFit &operator=(const RotationalFit &);

  public:
    // construct Reference for reference molecule, then
    // construct the RotationalFit.
    /**
     * RotationalFit constructor. It takes a reference to which a system 
     * can be fitted
     */
    RotationalFit(Reference *);
    /**
     * RotationalFit deconstructor
     */
    ~RotationalFit();
    /**
     * Method to fit your System to the Reference
     */
    void fit(gcore::System *)const;

    /**
     * accessor to the reference;
     */
    Reference * getReference() { return d_ref; }
    
    struct Exception: public gromos::Exception{
      Exception(const std::string &what): 
	gromos::Exception("RotationalFit",what){}
    };
    
  };

}

#endif    
