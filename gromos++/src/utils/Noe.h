#ifndef INCLUDED_UTILS_NOE
#define INCLUDED_UTILS_NOE

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gcore{
  class System;
}


namespace utils{
  class Noe_i;
  /**
   * Class Noe
   * The Noe class stores and analyses Noe information
   *
   * It stores the virtual atoms that define the NOE distance and the 
   * distances
   *
   * @class Noe
   * @author R. Buergi and M.A. Kastenholz
   * @ingroup utils
   * @todo finish the documentation
   */  
  class Noe{

    Noe_i *d_this;

    // not implemented
    Noe();
    Noe(const Noe &);
    Noe &operator=(const Noe&);
    
  
  public:
    Noe(const gcore::System &sys, const string &line);
  
    double distance(int i)const;
    // distance including correction.
    double correctedReference(int i)const;

    // Distance restraint for GROMOS96 compatibility mode... ;)
    string distRes(int i)const;

    // info string
    string info(int i)const;


    // accessors
    int numDistances()const;
    int numReferences()const;
    double reference(int i)const;
    double correction(int type);

    //method
    void setcorrection(int type, double correction);

    
    struct Exception: public gromos::Exception{
      Exception(const string &str): gromos::Exception("Noe", str){}
    };
    
  };

    

} /* namespace */

#endif
