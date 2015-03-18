// gio_OutBuildingBlock.h
#ifndef INCLUDED_OUTBUILDINGBLOCK
#define INCLUDED_OUTBUILDINGBLOCK

#include<string>
#include <set>

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
    /**
     * write a single building block 
     * @param bb the building block
     * @param type the type (either solute or end group)
     */
    void writeSingle(const gcore::BbSolute & bb, BBType type);
    /**
     * write a single solvent building block
     * @param bb the solvent building block
     */
    void writeSolvent(const gcore::SolventTopology & bb);
  };
}

#endif

