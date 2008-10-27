// TrajArray.h
#ifndef TRAJARRAY_H
#define TRAJARRAY_H
#include "../gcore/Box.h"
#include "../gcore/Molecule.h"
#include "../gcore/System.h"

namespace utils {
  /**
   * Class TrajArray
   * The TrajArray contains a complete trajectory
   *
   * If the coordinates of multiple frames is needed, the TrajArray can
   * be used to access these simultaneously without needing to store
   * the topological information of your system hundreds of times
   *
   * @class TrajArray
   * @author T. Hansson and V. Kraeutler
   * @ingroup utils
   */
  class TrajArray {
  public:
    /**
     * Constructor
     */
    TrajArray(const gcore::System &sys);

    /**
     * Destructor
     **/
    ~TrajArray();

    /**
     * store a frame in the array
     * @param sys the system to store
     * @param frameIndex index of the frame
     */
    void store(const gcore::System &sys,
        const unsigned int frameIndex);

    /**
     * extract a system from the array
     * @param sys the resulting system
     * @param frameIndex the index of the frame you want to extract
     */
    void extract(gcore::System &sys,
        const unsigned int frameIndex) const;

    /**
     * the number of of coordinates stored per frame
     */
    inline unsigned int numAtoms();

  protected:
    std::vector<double *> trajectoryData;
    // number of coordinates per frame
    unsigned int nAtoms;

  };
}
#endif                                                                 
