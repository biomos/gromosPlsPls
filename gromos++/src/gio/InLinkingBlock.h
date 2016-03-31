// gio_InBuildingBlock.h

#ifndef INCLUDED_GIO_IBUILDING
#define INCLUDED_GIO_IBUILDING

#include <string>
#include "../gromos/Exception.h"
#include "../gcore/BbLink.h"

namespace gcore{
  class BuildingBlock;
}

namespace gio{

  class InLinkingBlock_i;
  /**
   * Class InLinkingBlock
   * Instream that can read in a topology linking file
   *
   * @class InLinkingBlock
   * @author C. Oostenbrink
   * @ingroup gio
   * @sa gcore::BuildingBlock
   * @sa gcore::BbSolute
   + @sa gcore::BbLink
   * @sa gcore::SolventTopology
   */
  class InLinkingBlock{
    InLinkingBlock_i *d_this;
    // prevent default construction, copying and assignment
    InLinkingBlock();
    InLinkingBlock(const InLinkingBlock &);
    InLinkingBlock &operator=(const InLinkingBlock &);
    
  public:
    /**
     * open a building block file
     * @param file the file
     */
    InLinkingBlock(std::string file);
    
    ~InLinkingBlock();

    /**
     * access the building blocks
     */
    const std::vector<gcore::BbLink> &links()const;


    /**
     * access the force field code
     */
    const std::string forceField()const;

    /**
     * access the title
     */
    const std::string title()const;

    /**
     * @struct Exception
     * @ingroup gio::InBuildingBlock
     * The exception type for BB input exceptions
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what) : 
	gromos::Exception("InLinkingBlock", what){}
    };
  };
}
#endif
