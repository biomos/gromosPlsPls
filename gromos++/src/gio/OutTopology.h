// gio_OutTopology.h

#ifndef INCLUDED_STRING
#include<string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_SET
#include <set>
#define INCLUDED_SET
#endif

namespace gcore{
  class System;
  class GromosForceField;
}

namespace gio{
  /**
   * Class OutTopology
   * an outstream that defines how a gromos96 topology should be written out
   * 
   * The outputted topology does not have all the comments that a gromos96
   * topology has, but it can be read by GROMOS96 without problem
   *
   * @class OutTopology
   * @author B.C. Oostenbrink (actually, I am not sure)
   * @ingroup gio
   * @todo finish documentation
   */
  class OutTopology{
    std::string d_title;
    std::ostream &d_os;
    
    // not Implemented
    OutTopology();
    OutTopology(const OutTopology&);
    OutTopology &operator=(const OutTopology&);
  public:
    OutTopology(std::ostream &os);
    ~OutTopology();
    void setTitle(const std::string &title);
    void write(const gcore::System &sys, const gcore::GromosForceField &gff);
    
    //    OutTopology &operator<<(const gcore::Simulation &sim);
  };
}
