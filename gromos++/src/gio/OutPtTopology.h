// gio_OutTopology.h

#ifndef INCLUDED_OUTTOPOLOGY
#define INCLUDED_OUTTOPOLOGY "INCLUDED_OUTTOPOLOGY"
#ifndef INCLUDED_STRING
#include<string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_SET
#include <set>
#define INCLUDED_SET
#endif

namespace gcore{
  class PtTopology;
  class System;
}

namespace gio{
  /**
   * Class OutPtTopology
   * an out stream that defines how a perturbation topology should be written out
   *
   * @class OutPtTopology
   * @author N. Schmid
   * @ingroup gio
   */
  class OutPtTopology{
    std::string d_title;
    std::ostream &d_os;
  public:
    /**
     * construct the OutPtTopology from a stream
     * @param os the out stream to which the perturbation topology is written
     */
    OutPtTopology(std::ostream &os);
    /**
     * destruct the OutPtTopology
     */
    ~OutPtTopology();
    /**
     * sets the content of the title block
     * @param title the title
     */
    void setTitle(const std::string &title);
    /**
     * writes out a perturbation topology between two states A and B.
     * @param pt the perturbation topology
     * @param sys the reference system to detect the hydrogens. If NULL 
     *        (default) partitioning into hydrogen and non-hydrogen blocks
     *        is omitted.
     * @param a the index of state A (default 0)
     * @param b the index of state B (default 1)
     */
    void write(const gcore::PtTopology & pt,
               gcore::System * sys = NULL, int a=0, int b=1);
    /**
     * writes out a multiple perturbation topology (nonbonded only)
     * @param pt the perturbation topology
     */
    void write_multiple(const gcore::PtTopology & pt);
  };
}

#endif

