// gio_OutBuildingBlock.h

#ifndef INCLUDED_STRING
#include<string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_SET
#include <set>
#define INCLUDED_SET
#endif

#ifndef INCLUDED_OUTBUILDINGBLOCK
#define INCLUDED_OUTBUILDINGBLOCK

namespace gcore{
  class BuildingBlock;
  class BbSolute;
}

namespace gio{
  /**
   * Class OutBuildingBlock
   * an outstream that defines how a GROMOS molecular building block file should be written out
   *
   * @class OutBuildingBlock
   * @author N. Schmid
   * @ingroup gio
   */
  class OutBuildingBlock{
    std::string d_title;
    std::ostream &d_os;

    // prevent copying and assignment
    OutBuildingBlock();
    OutBuildingBlock(const OutBuildingBlock&);
    OutBuildingBlock &operator=(const OutBuildingBlock&);
  public:
    /**
     * The type of a building block
     */
    enum BBType { BBTypeSolute, BBTypeSolvent, BBTypeEnd };
    /**
     * create from a stream
     */
    OutBuildingBlock(std::ostream &os);
    ~OutBuildingBlock();
    /**
     * set the title
     * @param title the title
     */
    void setTitle(const std::string &title);
    /**
     * write the BuildingBlock as a molecular topology building block file
     * @param bb the building block data
     */
    void write(const gcore::BuildingBlock &bb);
  private:
    void writeSingle(const gcore::BbSolute & bb, BBType type);
    void writeSolvent(const gcore::SolventTopology & bb);
  };
}

#endif

