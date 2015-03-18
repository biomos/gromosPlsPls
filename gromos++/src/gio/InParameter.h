// gio_InParameter.h

#ifndef INCLUDED_GIO_IPARAMETER
#define INCLUDED_GIO_IPARAMETER

#include <string>
#include "../gromos/Exception.h"

namespace gcore{
  class GromosForceField;
}

namespace gio{

  class InParameter_i;
  /**
   * Class InParameter
   * defines an instream that can read a GROMOS ifp-file
   *
   * The GROMOS ifp file is read in and stored in a GromosForceField
   * format. This means that VdW-parameters are already calculated to
   * the individual pairs, taking all special tables in the manual into 
   * account.
   *
   * @class InParameter
   * @author B.C. Oostenbrink
   * @ingroup gio
   * @sa gcore::GromosForceField
   */
  class InParameter{
    InParameter_i *d_this;
    // prevent default construction, copying and assignment
    InParameter();
    InParameter(const InParameter &);
    InParameter &operator=(const InParameter &);
    
  public:
    /**
     * open a parameter file
     * @param file the file to open
     */
    InParameter(std::string file);
    
    ~InParameter();

    /**
     * access to the force field read
     */
    const gcore::GromosForceField &forceField()const;

    /**
     * access to the title
     */
    const std::string title()const;

    /**
     * The exception type for parameter reading
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what) : 
	gromos::Exception("InParameter", what){}
    };
  };
}
#endif
