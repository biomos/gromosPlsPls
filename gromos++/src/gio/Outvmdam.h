// gio_Outvmdam.h

#ifndef INCLUDED_GIO_OUTvmdam
#define INCLUDED_GIO_OUTvmdam

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
  class Outvmdam_i;
  
  class Outvmdam: public OutCoordinates{
    Outvmdam_i *d_this;
    // not implemented
    Outvmdam(const Outvmdam &);
    Outvmdam &operator=(const Outvmdam&);
  public:
    Outvmdam();
    Outvmdam(std::ostream &os);
    ~Outvmdam();
    void select(const std::string &thing);
    void open(std::ostream &os);
    void close();
    void writeTitle(const std::string &title);
    Outvmdam &operator<<(const gcore::System &sys);
  };
}

#endif
