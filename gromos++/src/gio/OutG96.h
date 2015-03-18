// gio_OutG96.h

#ifndef INCLUDED_GIO_OUTG96
#define INCLUDED_GIO_OUTG96

#include<string>
#include "OutCoordinates.h"

namespace gcore{
  class System;
  class Box;
}

namespace utils {
  class AtomSpecifier;
}

namespace gio{
  class OutG96_i;
  /**
   * Class OutG96
   * is of type OutCoordinates and defines how a gromos96 trajectory is 
   * printed out (POSITIONRED block) 
   *
   * @class OutG96
   * @author R. Buergi
   * @ingroup gio
   */
  class OutG96: public OutCoordinates{
    OutG96_i *d_this;
    // prevent copying and assignment
    OutG96(const OutG96 &);
    OutG96 &operator=(const OutG96&);
  public:
    OutG96();
    OutG96(std::ostream &os);
    ~OutG96();
    void select(const std::string &thing);
    void open(std::ostream &os);
    void close();
    void writeTitle(const std::string &title);
    void writeTimestep(const int step, const double time);
    void writeGenBox(const gcore::Box &box);
    void writeTriclinicBox(const gcore::Box &box);
    OutG96 &operator<<(const gcore::System &sys);
    OutG96 &operator<<(const utils::AtomSpecifier & atoms);
  };
}

#endif
