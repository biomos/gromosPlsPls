// gio_InG96.h

#ifndef INCLUDED_GIO_INCOORD
#define INCLUCED_GIO_INCOORD

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gcore{
  class System;
}

namespace gio{
  
  class InG96_i;
  /**
   * Class InG96
   * Defines an instream that can read any GROMOS96 coordinate file
   * or trajectory file.
   *
   * The instream can handle POSITION POSITIONRED VELOCITY and VELOCITYRED 
   * blocks in gromos96 files
   *
   * @class InG96
   * @author R. Buergi
   * @author M.K. Kastenholz, B.C. Oostenbrink (solvent)
   * @author M. Christen (velocities)
   * @sa gcore::System
   * @todo finish documentation
   */
  class InG96{
    InG96_i *d_this;
    // not implemented
    InG96(const InG96&);
    InG96 &operator=(const InG96&);
  public:
    // Constructors
    InG96();
    InG96(const string &name);
    ~InG96();

    // Methods
    void select(const string &thing);
    void open(const string &name);
    void close();
    InG96 &operator>>(gcore::System &sys);

    // Accessors
    const string &title()const;
    const string &name()const;
    bool eof()const;
    
    //Exceptions
    struct Exception: public gromos::Exception{
      Exception(const string& what_arg) : gromos::Exception("InG96", what_arg){}
    };
  };
}
#endif
