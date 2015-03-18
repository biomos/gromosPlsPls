// gio_Outvmdam.h

#ifndef INCLUDED_GIO_OUTvmdam
#define INCLUDED_GIO_OUTvmdam

#include<string>
#include "OutCoordinates.h"

namespace gcore{
  class System;
}

namespace gio{
  class Outvmdam_i;
  /**
   * Class Outvmdam
   * is of type OutCoordinates and defines how a trajectory should be
   * written out in "Amber Coordinates" that can be read by VMD.
   *
   * @class Outvmdam
   * @author M.K. Kastenholz
   * @ingroup gio
   */
  class Outvmdam: public OutCoordinates{
    Outvmdam_i *d_this;
    double factor;
    // prevent copying and assignment
    Outvmdam(const Outvmdam &);
    Outvmdam &operator=(const Outvmdam&);
  public:
    Outvmdam(double factor = 10.0);
    Outvmdam(std::ostream &os, double factor = 10.0);
    ~Outvmdam();
    void select(const std::string &thing);
    void open(std::ostream &os);
    void close();
    void writeTitle(const std::string &title);
    void writeTimestep(const int step, const double time);
    Outvmdam &operator<<(const gcore::System & sys);
    Outvmdam &operator<<(const utils::AtomSpecifier & atoms);
  };
}

#endif
