// gio_OutTopology.h

#ifndef INCLUDED_STRING
#include<string>
#define INCLUDED_STRING
#endif


namespace gcore{
  class System;
  class GromosForceField;
}

namespace gio{
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
