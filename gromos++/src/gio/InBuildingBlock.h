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
  /**
   * Class InBuildingBlock
   * Instream that can read in a gromos96 mtb-file
   *
   * In addition to the standard blocks in the mtb-file, this stream can also
   * read in MTBUILDBLEND blocks which contain the definition of end-groups
   * These blocks should appear at the end of the file, after the 
   * MTBUILDBLSOLVENT blocks
   *
   * @class InBuildingBlock
   * @author B.C. Oostenbrink
   * @ingroup gio
   * @sa gcore::BuildingBlock
   * @sa gcore::BbSolute
   * @sa gcore::SolventTopology
   * @todo finish documentation
   */
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
