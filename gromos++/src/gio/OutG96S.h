// gio_OutG96S.h

#ifndef INCLUDED_GIO_OUTG96S
#define INCLUDED_GIO_OUTG96S

#include<string>
#include "OutCoordinates.h"

namespace gcore{
  class System;
}

namespace utils{
  class AtomSpecifier;
}

namespace gio{
  class OutG96S_i;
  /**
   * Class OutG96S
   * is of type OutCoordinates and defines how a single coordinate file 
   * is printed out (POSITION block)
   *
   * @class OutG96S
   * @author R. Buergi
   * @author M.K. Kastenholz, B.C. Oostenbrink (solvent)
   * @ingroup gio
   */
  class OutG96S: public OutCoordinates{
    OutG96S_i *d_this;
    bool posres;
    // prevent copying and assignment
    OutG96S(const OutG96S &);
    OutG96S &operator=(const OutG96S&);
  public:
    OutG96S(bool posres = false);
    OutG96S(std::ostream &os, bool posres = false);
    ~OutG96S();
    void select(const std::string &thing);
    void open(std::ostream &os);
    void close();
    void writeTitle(const std::string &title);
    void writeTimestep(const int step, const double time);
    OutG96S &operator<<(const gcore::System & sys);
    OutG96S &operator<<(const utils::AtomSpecifier & atoms);
  };
}

#endif
