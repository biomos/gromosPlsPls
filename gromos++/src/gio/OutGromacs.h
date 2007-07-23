// gio_OutGromacs.h

#ifndef INCLUDED_STRING
#include<string>
#define INCLUDED_STRING
#endif


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
   * @todo finish documentation
   */
  class OutGromacs{
    std::string d_title;
    std::ostream &d_os;
    // not Implemented
    OutGromacs();
    OutGromacs(const OutGromacs&);
    OutGromacs &operator=(const OutGromacs&);
  public:
    OutGromacs(std::ostream &os);
    ~OutGromacs();
    void setTitle(const std::string &title);
    void write(const gcore::System &sys, const gcore::GromosForceField &gff);
    
    //    OutTopology &operator<<(const gcore::Simulation &sim);
  };
}
