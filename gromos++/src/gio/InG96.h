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
