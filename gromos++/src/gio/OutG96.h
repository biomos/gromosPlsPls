// gio_OutG96.h

#ifndef INCLUDED_GIO_OUTG96
#define INCLUDED_GIO_OUTG96

#ifndef INCLUDED_STRING
#include<string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_GIO_OUTCOORDINATES
#include "OutCoordinates.h"
#endif

namespace gcore{
  class System;
}

namespace gio{
  class OutG96_i;
  
  class OutG96: public OutCoordinates{
    OutG96_i *d_this;
    // not implemented
    OutG96(const OutG96 &);
    OutG96 &operator=(const OutG96&);
  public:
    OutG96();
    OutG96(std::ostream &os);
    ~OutG96();
    void select(const std::string &thing);
    void open(std::ostream &os);
    void close();
    void writeTitle(const std::string &title);
    OutG96 &operator<<(const gcore::System &sys);
  };
}

#endif
