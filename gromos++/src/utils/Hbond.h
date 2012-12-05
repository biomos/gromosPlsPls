/** 
 * @file:   Hbond.h
 * Author: siggj
 *
 * Created on October 10, 2012, 10:15 AM
 */

#ifndef INCLUDED_HBOND_CALC
#define	INCLUDED_HBOND_CALC

#include "Hbond_calc_2c.h"
#include "Hbond_calc_3c.h"


namespace utils {

  /**
   * Class HB
   * purpose: To serve as an interface between the hbond.cc and the 2- or 3-centred H-Bond-
   * calculations. 
   * @class HB
   * @author J.Sigg
   * @ingroup utils
   */
  class HB {
    HB2c_calc hb2c_calc;
    HB3c_calc hb3c_calc;
    HBPara2c hbpara2c;
    HBPara3c hbpara3c;
    bool do3c, do_native;
  public:
    /**
     * Method to initialize the calculation.
     */
    void init();
    /**
     * Method, that chooses, which calculation have to be done.
     */
    void calc();
    /**
     * Method to clear the H-bonds after calculating the native H-bonds.
     */
    void clear();
    /**
     * Method that stores the time of the trajectories to the H-bonds.
     * @param times
     */
    void settime(double times);
    /**
     * Method that prints out all information about the H-bonds.     
     */
    void printstatistics();
    /**
     * Constructor which gives the parameters to the belonging class.
     * @param sys
     * @param args
     * @param hbparas2c
     * @param hbparas3c
     */
    HB(gcore::System &sys, args::Arguments &args, HBPara2c hbparas2c, HBPara3c hbparas3c);

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

