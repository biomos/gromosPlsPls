// gio_InTopology.h

#ifndef INCLUDED_GIO_ITOPOLOGY
#define INCLUDED_GIO_ITOPOLOGY

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gcore{
  class System;
  class GromosForceField;
}

namespace gio{

  class InTopology_i;
  /**
   * Class InTopology
   * defines an instream that can read in a GROMOS96 topology
   *
   * The data that is read in is split up into a System and a
   * GromosForceField
   *
   * @class InTopology
   * @author R. Buergi
   * @author B.C. Oostenbrink (massType, Solvent)
   * @sa gcore::System
   * @sa gcore::GromosForceField
   * @todo finish documentation
   */
  class InTopology{
    InTopology_i *d_this;
    // not implemented
    InTopology();
    InTopology(const InTopology &);
    InTopology &operator=(const InTopology &);
    
  public:
    // Constructors
    InTopology(string str);
    
    ~InTopology();

    // methods
    const gcore::System &system()const;
    const gcore::GromosForceField &forceField()const;

    // accessors
    const string &version()const;
    const string &title()const;

    //Exceptions
    struct Exception: public gromos::Exception{
      Exception(const string& what) : gromos::Exception("InTopology", what){}
    };
  };
}
#endif
