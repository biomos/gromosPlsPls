// gio_OutG96S.h

#ifndef INCLUDED_GIO_OUTG96S
#define INCLUDED_GIO_OUTG96S

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
  class OutG96S_i;
  
  class OutG96S: public OutCoordinates{
    OutG96S_i *d_this;
    // not implemented
    OutG96S(const OutG96S &);
    OutG96S &operator=(const OutG96S&);
  public:
    OutG96S();
    OutG96S(std::ostream &os);
    ~OutG96S();
    void select(const std::string &thing);
    void open(std::ostream &os);
    void close();
    void writeTitle(const std::string &title);
    OutG96S &operator<<(const gcore::System &sys);
  };
}

#endif
