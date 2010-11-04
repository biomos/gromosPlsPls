// gio_OutPdb.h

#ifndef INCLUDED_GIO_OUTPdb
#define INCLUDED_GIO_OUTPdb

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
  class OutPdb_i;
  /**
   * Class OutPdb
   * is of type OutCoordinates and defines how a PDB-file is written out
   * 
   * @class OutPdb
   * @author R. Buergi
   * @author M.K. Kastenholz, B.C. Oostenbrink (solvent)
   * @ingroup gio
   */
  class OutPdb: public OutCoordinates{
    OutPdb_i *d_this;
    // prevent copying and assignment
    OutPdb(const OutPdb &);
    OutPdb &operator=(const OutPdb&);
    double factor;
  public:
    OutPdb(double factor = 10.0);
    OutPdb(std::ostream &os, double factor = 10.0);
    ~OutPdb();
    void select(const std::string &thing);
    void open(std::ostream &os);
    void close();
    void writeTitle(const std::string &title);
    void writeTimestep(const int step, const double time);
    OutPdb &operator<<(const gcore::System &sys);
    OutPdb &operator<<(const utils::AtomSpecifier & atoms);
  };
}

#endif
