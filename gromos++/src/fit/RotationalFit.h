// fit_RotationalFit.h

#ifndef INCLUDED_FIT_ROTATIONALFIT
#define INCLUDED_FIT_ROTATIONALFIT


#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gcore{
  class System;
}

namespace fit{
  class RotationalFit_i;
  class Reference;

  class RotationalFit{
    const Reference *d_ref;

    // not implemented
    RotationalFit();
    RotationalFit(const RotationalFit &);
    RotationalFit &operator=(const RotationalFit &);

  public:
    // construct Reference for reference molecule, then
    // construct the RotationalFit.
    RotationalFit(Reference *);
    ~RotationalFit();
    
    void fit(gcore::System *)const;
    
    struct Exception: public gromos::Exception{
      Exception(const string &what): gromos::Exception("RotationalFit",what){}
    };
    
  };

}

#endif    
