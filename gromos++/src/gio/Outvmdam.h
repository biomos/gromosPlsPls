// gio_Outvmdam.h

#ifndef INCLUDED_GIO_OUTvmdam
#define INCLUDED_GIO_OUTvmdam

#ifndef INCLUDED_STRING
#include<string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_GIO_OUTCOORDINATES
#include "OutCoordinates.h"
#endif

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
    // prevent copying and assignment
    Outvmdam(const Outvmdam &);
    Outvmdam &operator=(const Outvmdam&);
  public:
    Outvmdam();
    Outvmdam(std::ostream &os);
    ~Outvmdam();
    void select(const std::string &thing);
    void open(std::ostream &os);
    void close();
    void writeTitle(const std::string &title);
    void writeTimestep(const int step, const double time);
    Outvmdam &operator<<(const gcore::System &sys);
  };
}

#endif
