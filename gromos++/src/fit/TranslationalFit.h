// fit_TranslationalFit.h

#ifndef INCLUDED_FIT_TRANSLATIONALFIT
#define INCLUDED_FIT_TRANSLATIONALFIT


namespace gcore{
  class System;
}

namespace gmath{
  class Vec;
}

namespace fit{
  class TranslationalFit_i;
  class Reference;

  class TranslationalFit{
    TranslationalFit_i *d_this;

    // not implemented
    TranslationalFit();
    TranslationalFit(const TranslationalFit &);
    TranslationalFit &operator=(const TranslationalFit &);

  public:
    // construct Reference for reference molecule, then
    // construct the TranslationalFit.
    TranslationalFit(Reference *);
    ~TranslationalFit();
    
    //    void fitToCom(gcore::System *)const;
    void fit(gcore::System *)const;

    const gmath::Vec &com()const;
    const gmath::Vec &cog()const;
    
  };

}

#endif    
