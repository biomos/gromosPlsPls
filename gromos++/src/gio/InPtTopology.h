// gio_InTopology.h

#ifndef INCLUDED_GIO_INPTTOPOLOGY
#define INCLUDED_GIO_INPTTOPOLOGY

#include <string>
#include "../gromos/Exception.h"

namespace gcore{
  class System;
  class PtTopology;
}

namespace gio{
  class InPtTopology_i;
  /**
   * Class InPtTopology
   * defines an instream that can read in a perturbation topology
   *
   * @class InPtTopology
   * @ingroup gio
   * @author B.C. Oostenbrink  N. Schmid
   * @sa gcore::System
   * @sa gcore::PtTopology
   */
  class InPtTopology{
    InPtTopology_i *d_this;
    
  public:
    /**
     * Construct the InPtTopology from a file name and parse the content
     */
    InPtTopology(std::string str);
    /**
     * Destructor
     */
    ~InPtTopology();
    
    /**
     * Accessor to the perturbation topology read
     */
    const gcore::PtTopology &ptTopo()const;
    /**
     * Accessor to the title
     */
    const std::string title()const;

    /**
     * @struct Exception
     * The Exception that occurs when reading a perturbation topology
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what) : 
	gromos::Exception("InPtTopology", what){}
    };
  };
}
#endif
