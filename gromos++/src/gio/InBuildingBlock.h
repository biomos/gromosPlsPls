// gio_InBuildingBlock.h

#ifndef INCLUDED_GIO_IBUILDING
#define INCLUDED_GIO_IBUILDING

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gcore{
  class BuildingBlock;
}

namespace gio{

  class InBuildingBlock_i;

  class InBuildingBlock{
    InBuildingBlock_i *d_this;
    // not implemented
    InBuildingBlock();
    InBuildingBlock(const InBuildingBlock &);
    InBuildingBlock &operator=(const InBuildingBlock &);
    
  public:
    // Constructors
    InBuildingBlock(string str);
    
    ~InBuildingBlock();

    // methods
    const gcore::BuildingBlock &building()const;

    // accessors
    const string &title()const;

    //Exceptions
    struct Exception: public gromos::Exception{
      Exception(const string& what) : gromos::Exception("InBuildingBlock", what){}
    };
  };
}
#endif
