// gio_InTopology.h

#ifndef INCLUDED_GIO_INTOPOLOGY
#define INCLUDED_GIO_INTOPOLOGY

#include <string>
#include "../gromos/Exception.h"

namespace gcore{
  class System;
  class GromosForceField;
}

namespace gio{
  class InTopology_i;
  /**
   * Class InTopology
   * defines an instream that can read in a GROMOS topology
   *
   * The data that is read in is split up into a System and a
   * GromosForceField
   *
   * @class InTopology
   * @ingroup gio
   * @author R. Buergi
   * @author B.C. Oostenbrink (massType, Solvent)
   * @sa gcore::System
   * @sa gcore::GromosForceField
   */
  class InTopology{
    InTopology_i *d_this;
    
  public:
    /**
     * open a topology file
     * @param file the topology file to open
     */
    InTopology(std::string file);
    ~InTopology();
    
    /**
     * access to the system that was read
     */
    const gcore::System &system()const;
    /**
     * access to the force field that was read
     */
    const gcore::GromosForceField &forceField()const;

    /**
     * access to the version string
     */
    const std::string &version()const;
    /**
     * access to the title
     */
    const std::string title()const;

    /**
     * The exception type for toplogy reading
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what) : 
	gromos::Exception("InTopology", what){}
    };
  };
}
#endif
