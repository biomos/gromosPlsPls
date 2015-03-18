// gio_OutGromacs.h

#ifndef INCLUDED_OUTGROMACS
#define INCLUDED_OUTGROMACS

#include<string>


namespace gcore{
  class System;
  class GromosForceField;
}

namespace gio{
  /**
   * Class OutGromacs
   * an outstream that writes a topology that looks similar to a 
   * gromacs topology
   * 
   * @class OutGromacs
   * @author B.C. Oostenbrink
   * @ingroup gio
   */
  class OutGromacs{
    std::string d_title;
    std::ostream &d_os;
    // prevent copying and assignment
    OutGromacs();
    OutGromacs(const OutGromacs&);
    OutGromacs &operator=(const OutGromacs&);
  public:
    /**
     * construct using an output stream
     */
    OutGromacs(std::ostream &os);
    ~OutGromacs();
    /**
     * set the title
     * @param title the title
     */
    void setTitle(const std::string &title);
    /**
     * write the system and force-field parameters in gromacs format
     * @param sys the system
     * @param ggf the force field
     */
    void write(const gcore::System &sys, const gcore::GromosForceField &gff);
  };
}

#endif
