// gio_InParameter.h

#ifndef INCLUDED_GIO_IPARAMETER
#define INCLUDED_GIO_IPARAMETER

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gcore{
  //  class System;
  class GromosForceField;
}

namespace gio{

  class InParameter_i;
  /**
   * Class InParameter
   * defines an instream that can read a GROMOS96 ifp-file
   *
   * The GROMOS96 ifp file is read in and stored in a GromosForceField
   * format. This means that vdw-parameters are already calculated to 
   * the individual pairs, taking all special tables in the manual into 
   * account'
   *
   * @class InParameter
   * @author B.C. Oostenbrink
   * @ingroup gio
   * @sa gcore::GromosForceField
   * @todo finish documentation
   */
  class InParameter{
    InParameter_i *d_this;
    // not implemented
    InParameter();
    InParameter(const InParameter &);
    InParameter &operator=(const InParameter &);
    
  public:
    // Constructors
    InParameter(string str);
    
    ~InParameter();

    // methods
    //const gcore::System &system()const;
    const gcore::GromosForceField &forceField()const;

    // accessors
    //const string &version()const;
    const string &title()const;

    //Exceptions
    struct Exception: public gromos::Exception{
      Exception(const string& what) : gromos::Exception("InParameter", what){}
    };
  };
}
#endif
