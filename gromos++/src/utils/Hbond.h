/**
 * @file:   Hbond.h
 * Author: siggj
 *
 * Created on October 10, 2012, 10:15 AM
 */

#ifndef INCLUDED_HBOND_CALC
#define	INCLUDED_HBOND_CALC

#include <vector>

#ifdef OMP
#include <omp.h>
#endif

#include "AtomSpecifier.h"

#include "Hbond_calc_2c.h"
#include "Hbond_calc_3c.h"
#include "Hbond_calc_bridges.h"
#include "CubeSystem.hcc"


namespace utils {

  /**
   * Class HB
   * purpose: To serve as an interface between the hbond.cc and the 2-, 3-centred, or solvent bridges H-Bond-
   * calculations.
   * @author J.Sigg, M.Setz
   * @ingroup utils
   * @class HB
   */
  class HB {
    bool do3c, do_native, sort_occ, doBridges, reduce, vacuum;
    double higherthan;
    HB2c_calc hb2c_calc;
    HB3c_calc hb3c_calc;
    HB_bridges hb_bridges;
    HBPara2c hbpara2c;
    HBPara3c hbpara3c;

  public:

    /**
     * Method that calls the requested calculations.
     */
    void calc(CubeSystem<int>&, CubeSystem<int>&,  CubeSystem<Key2c>&);
    /**
     * Method to clear the H-bonds after calculating the native H-bonds.
     */
    void clear();
    /**
     * Method that stores the time of the trajectories to the H-bonds.
     */
    void settime(double times);
    /**
     * Method that prints out all information about the H-bonds.
     */
    void printstatistics();
    /**
     * Constructor which gives the parameters to the belonging class.
     */
    HB(gcore::System &sys, args::Arguments &args, HBPara2c hbparas2c, HBPara3c hbparas3c);

    //function to merge hbond maps into one output vector: for omp parallelized trajectories
    /**
     * Method that merges all H-bond objects for OpenMP parallelized trajectories.
     */
    void merge(HB&);

    /**
     * Method that prepares everything for native H-bond calculation.
     */
    void prepare_native(CubeSystem<int>&, CubeSystem<int>&, CubeSystem<Key2c>&);

    static HBPara2c mk_hb2c_paras(const vector<double> &hbparas);
    static HBPara3c mk_hb3c_paras(const vector<double> &hbparas);

    enum arg_name {
      DISTANCE = 0,
      ANGLE = 1,
      SUM = 2,
      DIHEDRAL = 3
    };
  }; //end class HB




}
#endif	/* NEWFILE_H */

