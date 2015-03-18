/**
 * @file InBFactorOccupancy.h
 * a file to read B-factor and occupancy data
 */

#ifndef INCLUDED_INBFACTOROCCUPANCY_H
#define	INCLUDED_INBFACTOROCCUPANCY_H

#include <string>
#include <vector>

#include "Ginstream.h"

namespace gio {
  /**
   * @struct BFactorOccupancyData
   * @ingroup gio
   * @author N.Schmid, F. Freitag
   *
   * A class to hold B-factor and occupancy data.
   */
  struct BFactorOccupancyData {
    /**
     * the B factor
     */
    double b_factor;
    /**
     * the occupancy value
     */
    double occupancy;
  };

  /**
   * @class InBFactorOccupancy
   * @ingroup gio
   * @author N. Schmid, F. Freitag
   * @brief reads a B-factor and Occupancy file
   * A class to read files containing B-factor and occupancy information
   *
   * Format of the B-factor and occupancy file:
   * The B-factors have to be given @f$\mathrm{nm}^2@f$.
   * @verbatim
TITLE
B-factors and occupancies for all atoms
END
BFACTOROCCUPANCY
# B-factor Occupancy
0.01  1.0
0.02  0.8
END
   @endverbatim
   */
  class InBFactorOccupancy {
  public:
    /**
     * constructor
     */
    InBFactorOccupancy() {}
    /**
     * construct from file name
     * @param file file name
     */
    InBFactorOccupancy(std::string file);
    /**
     * destructor
     */
    ~InBFactorOccupancy();
    /**
     * open a file
     * @param file file name
     */
    void open(std::string file);
    /**
     * close the file
     */
    void close();
    /**
     * get the mapping data
     * @return a map containing the IAC to element mapping
     */
    std::vector<BFactorOccupancyData> getData();
  protected:
    Ginstream file_stream;
  };
}

#endif	/* INCLUDED_INBFACTOROCCUPANCY_H */

