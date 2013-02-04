/**
 * @file:   Hbond_calc_2c.h
 * Author: siggj
 *
 * Created on October 10, 2012, 2:16 PM
 */

#ifndef INCLUDED_2C_H
#define	INCLUDED_2C_H

#include <fstream>
#include "Hbond_calc.h"


namespace utils {

  /**
   * Struct HBPara2c
   * purpose: to store all the paramaters (maximal distance, minimal angle) given
   * by the input file.
   * @struct HBPara2c
   * @author J.Sigg
   * @ingroup utils
   */
  struct HBPara2c {
    double maxdist, maxdist2, minangle;
  }; //end struct HBPara2c

  /**
   * Class HB2c
   * purpose: to store all 2-centred H-bonds, which have been found, with all 
   * the informations (distances, angles, indices) needed.
   * @class HB2c
   * @author J.Sigg
   * @ingroup utils
   */
  class HB2c {
    int num, index;
    double mean_dist, mean_angle;
    double dist, angle;
  public:
    /**
     * Method to clear, after the native calculation, all parameters already 
     * calculated.
     */
    void clear();
    /**
     * Method to calculate the mean distance and angle of an H-bond to print it out.
     */
    void calcmean();
    /**
     * Method which increases the number of the given H-bond by 1.
     */
    void incr_num();
    /**
     * Method to set the distance of the H-bond.
     * @param distance
     */
    void setdist(double &distance);
    /**
     * Method to set the angle of the H-bond.
     * @param ang
     */
    void setangle(double &ang);
    /**
     * Method, which gives back the number of times, the H-bond have been found.
     * @return 
     */
    int getnum();
    /**
     * Method, which gives the index of the H-bond back.
     * @return 
     */
    int getindex();
    /**
     * Method to set the index of the H-bond.
     * @param index
     */
    void setindex(int index);
    /**
     * Method, which gives the mean-distance back.
     * @return 
     */
    double getmeandist();
    /**
     * Method, which gives the mean-distance back.
     * @return 
     */
    double getmeanangle();

    /**
     * Constructor
     */
    HB2c() {
      dist = 0;
      angle = 0;
      num = 0;
    };
  }; //end class HB2c

  typedef std::map<unsigned int, HB2c> Hb2c_container;

  /**
   * Class HB2c_calc
   * purpose: Class, which inherit from HB_calc, to calculate the 2-centred H-bonds.
   * If a bond is found it is stored in the map of the class HB2c.
   * @class HB2c_calc
   * @author J. Sigg
   * @ingroup utils
   */
  class HB2c_calc : public HB_calc {
    Hb2c_container hb2cc;
    HBPara2c hbpara;
    std::vector<int> tsnum, tsnumHB;
    int numHb, maxindex;
    std::ofstream timeseriesHB, timeseriesHBtot;
    std::vector<double> tstime, tstimeHB;
    /**
     * Method that opens both timeseries files.
     * @param fi1
     * @param fi2
     */
    void opents(std::string fi1, std::string fi2);
    /**
     * Method to calculate if there is a H-bond between atom i and j.
     * @param i
     * @param j
     * @param bound_i
     */
    void calc2c(int i, int j, gmath::Vec &bound_i);

    /**
     * Method which gives the index from the donor back.
     * @param index
     * @return 
     */
    inline int getdon_ind(int index) {
      return index / acceptors.size();
    }

    /**
     * Method which gives the index from the acceptor back.
     * @param index
     * @return 
     */
    inline int getacc_ind(int index) {
      return index % acceptors.size();
    }
  public:
    /**
     * Constructor, which stores all parameters given from the input file.
     * @param para
     */
    HB2c_calc(HBPara2c para);

    /**
     * Destructor, which closes the timeseries files.
     */
    ~HB2c_calc() {
      timeseriesHB.close();
      timeseriesHBtot.close();
    }
    /**
     * Method to store the system file and the argument file for further use, and 
     * opens the timeseries file
     * @param sys
     * @param args
     */
    void setval(gcore::System &sys, args::Arguments &args);
    /**
     * Method to initialize the calculation.
     */
    void init();
    /**
     * Method to clear the parameters calculated during the native H-bond calculation.
     */
    void clear();
    /**
     * Method to calculate only H-bonds, which have been found during the native 
     * H-bond calculation.
     */
    void calc_native();
    /**
     * Method that loops over all atoms from DonorA and DonorB.
     */
    void calc2();
    /**
     * Method that prints out all 2-centered H-bonds.
     */
    void printstatistics();
    /**
     * Gives the occurrence of the H-Bond back.
     * @param index
     * @return 
     */
    int getnum(int index);
  }; //end class HB2c_calc

}
#endif	/* NEW_2C_H */
