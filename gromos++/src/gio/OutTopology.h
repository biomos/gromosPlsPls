// gio_OutTopology.h

#ifndef INCLUDED_OUTTOPOLOGY
#define INCLUDED_OUTTOPOLOGY

#include<string>
#include <set>

namespace gcore{
  class System;
  class GromosForceField;
}

namespace gio{
  /**
   * Class OutTopology
   * an outstream that defines how a GROMOS topology should be written out
   *
   * @class OutTopology
   * @author B.C. Oostenbrink
   * @ingroup gio
   */
  class OutTopology{
    std::string d_title;
    std::ostream &d_os;
    
    // prevent copying and assignment
    OutTopology();
    OutTopology(const OutTopology&);
    OutTopology &operator=(const OutTopology&);
  public:
    /**
     * create from a stream
     */
    OutTopology(std::ostream &os);
    ~OutTopology();
    /**
     * set the title
     * @param title the title
     */
    void setTitle(const std::string &title);
    /**
     * write the system and force field as a GROMOS topology to the stream
     * @param sys the system
     * @param ggf the force field
     */
    void write(const gcore::System &sys, const gcore::GromosForceField &gff);
    /**
     * write the system and force field as a GROMOS96 topology to the stream
     * @param sys the system
     * @param ggf the force field
     */
    void write96(const gcore::System &sys, const gcore::GromosForceField &gff);
  };
}
#endif
