// gio_InTopology.h

#ifndef INCLUDED_GIO_INPTTOPOLOGY
#define INCLUDED_GIO_INPTTOPOLOGY

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#define INCLUDED_GROMOS_EXCEPTION
#endif

namespace gcore{
  class System;
}

namespace gio{
  class InPtTopology_i;
  /**
   * Class InTopology
   * defines an instream that can read in a perturbation topology
   *
   * @class InTopology
   * @author B.C. Oostenbrink 
   * @sa gcore::System
   * @sa gcore::PtTopology
   * @todo finish documentation
   */
  class InPtTopology{
    InPtTopology_i *d_this;
    
  public:
    // Constructors
    InPtTopology(std::string str, int start=0);
    ~InPtTopology();
    
    // methods
    const gcore::PtTopology &ptTopo()const;

    // accessors
    const std::string title()const;

    //Exceptions
    struct Exception: public gromos::Exception{
      Exception(const std::string& what) : 
	gromos::Exception("InPtTopology", what){}
    };
  };
}
#endif
