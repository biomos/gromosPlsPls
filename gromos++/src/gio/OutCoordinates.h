// gio_OutCoordinates.h

#ifndef INCLUDED_GIO_OUTCOORDINATES
#define INCLUDED_GIO_OUTCOORDINATES

#ifndef INCLUDED_STRING
#include<string>
#define INCLUDED_STRING
#endif

namespace gcore{
  class System;
}

namespace gio{
  
  class OutCoordinates{
    // not implemented
    OutCoordinates(const OutCoordinates &);
    OutCoordinates &operator=(const OutCoordinates&);
  public:
    OutCoordinates(){}
    virtual ~OutCoordinates();
    virtual void open(std::ostream &os)=0;
    virtual void close()=0;
    virtual void writeTitle(const std::string &title)=0;
    virtual void select(const std::string &thing)=0;
    virtual OutCoordinates &operator<<(const gcore::System &sys)=0;
  };
}

#endif
