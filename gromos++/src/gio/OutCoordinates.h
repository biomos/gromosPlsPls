// gio_OutCoordinates.h

#ifndef INCLUDED_GIO_OUTCOORDINATES
#define INCLUDED_GIO_OUTCOORDINATES

#include<string>

namespace gcore {
  class System;
}

namespace utils {
  class AtomSpecifier;
}

namespace gio{
  /**
   * Class OutCoordinates
   * defines some basic features for an output stream that writes out
   * GROMOS coordinate or trajectory files
   *
   * @class OutCoordinates
   * @author R. Buergi
   * @ingroup gio
   */
  class OutCoordinates{
    // prevent copying and assignment
    OutCoordinates(const OutCoordinates &);
    OutCoordinates &operator=(const OutCoordinates&);
  public:
    OutCoordinates(){}
    virtual ~OutCoordinates();
    /**
     * open an output stream
     */
    virtual void open(std::ostream &os)=0;
    /**
     * close the output stream
     */
    virtual void close()=0;
    /**
     * write the title string
     */
    virtual void writeTitle(const std::string &title)=0;
    /**
     * write the time and step information
     */
    virtual void writeTimestep(const int step, const double time)=0;
    /**
     * select only parts of the system
     * @param thing ALL for everything, SOLVENT for solvent only, or anything else for solute
     */
    virtual void select(const std::string &thing)=0;
    /**
     * write a system to the stream
     */
    virtual OutCoordinates &operator<<(const gcore::System & sys)=0;
    /**
     * write an AtomSpecifier to the stream
     */
    virtual OutCoordinates &operator<<(const utils::AtomSpecifier & atoms)=0;
  };
}

#endif
