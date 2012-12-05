/** 
 * @file:   Hbond_calc_3c.h
 * Author: siggj
 *
 * Created on October 10, 2012, 2:17 PM
 */

#ifndef INCLUDED_3C_H
#define	INCLUDED_3C_H

#include <fstream>
#include "Hbond_calc.h"
#include "Hbond_calc_2c.h"


namespace utils {

  /**
   * Struct HBPara3c
   * purpose: to store all the paramaters (maximal distance, minimal angle, minimal
   * angle sum and maximal dihedral angle) given by the input file.
   * @struct HBPara3c
   * @author J.Sigg
   * @ingroup utils
   */
  struct HBPara3c {
    double maxdist, maxdist2, minangle;
    double minanglesum, maxdihedral;
  }; //end struct HBPara3c

  /**
   * Class HB3c
   * purpose: to store all 3-centred H-bonds, which have been found, with all 
   * the informations, (distances, angles, sum of angels, dihedral angles, 
   * indices) needed.
   * @class HB3c
   * @author J.Sigg
   * @ingroup utils
   */
  class HB3c {
    int num;
    double dihed_mean, angletot_mean;
    vector<double> mean_distance, mean_angle;
    vector<double> distance, angle;
    double dihed, angletot;
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
     * @param dist1
     * @param dist2
     */
    void setdist(double &dist1, double &dist2);
    /**
     * Method to set the angle and the sum of the angles of the H-bond.
     * @param angl1
     * @param angl2
     * @param angltot
     */
    void setangle(double &angl1, double &angl2, double &angltot);
    /**
     * Method to set the dihedral angle of the H-bond.
     * @param dihedral
     */
    void setdihed(double dihedral);
    /**
     * Method, which gives back the number of times, the H-bond have been found.
     * @return 
     */
    int getnum();
    /**
     * Method, which gives the mean-distance (first or second) back.
     * @param i
     * @return 
     */
    double getmeandist(int i);
    /**
     * Method, which gives the mean-angle (first or second) back.
     * @param i
     * @return 
     */
    double getmeanangle(int i);
    /**
     * Method, which gives the mean-sum of the angles back.
     * @return 
     */
    double getmeanangle_sum();
    /**
     * Method, which gives the mean-dihedral angle back.
     * @return 
     */
    double getmeandihedral();

    /**
     * Constructor
     */
    HB3c() {
      distance.resize(2);
      angle.resize(2);
      num = 0;
      dihed = 0.0;
      angletot = 0.0;
    };
  }; //end class HB3c

  typedef std::map<unsigned int, HB3c> Hb3c_container;

  /**
   * Class HB3c_calc
   * purpose: Class, which inherit from HB_calc, to calculate the 3-centred H-bonds.
   * If a bond is found it is stored in the map of the class HB3c.
   * @class HB3c_calc
   * @author J. Sigg
   * @ingroup utils
   */
  class HB3c_calc : public HB_calc {
    double min_angle_sum;
    double max_dihedral;
    std::ofstream timeseriesHB3c, timeseriesHB3ctot;
    std::vector<int> tsnum3c, tsnumHB3c;
    std::vector<double> tstime3c, tstimeHB3c;
    Hb3c_container hb3cc;
    HBPara3c hbpara;
    int numHb3c, index;
    /**
     * Method that opens both timeseries files.
     * @param fi1
     * @param fi2
     */
    void opents3c(std::string fi1, std::string fi2);
    /**
     * Method that loops over all k and calculate if there is a H-bond between 
     * atom i, j and k.
     * @param i
     * @param j
     * @param k_lim
     * @param bound_i
     */
    void calc3c(int i, int j, int k_lim, gmath::Vec &bound_i);

    /**
     * Method which gives the index from the donor back.
     * @param index
     * @return 
     */
    inline int getdon_ind(int index) {
      return (index / acceptors.size()) / acceptors.size();
    }

    /**
     * Method which gives the index from the first acceptor back.
     * @param index
     * @return 
     */
    inline int getacc_indj(int index) {
      return (index / acceptors.size()) % acceptors.size();
    }

    /**
     * Method which gives the index from the second acceptor back.
     * @param index
     * @return 
     */
    inline int getacc_indk(int index) {
      return index % acceptors.size();
    }

    /**
     * Method that gives back true if the sum of the angles between the atoms 
     * i, j and k is larger than the minimal sum of angles given.
     * @param i
     * @param angle_sum
     * @param acceptor_j
     * @param acceptor_k
     * @return 
     */
    inline bool anglesum(int i, double &angle_sum, gmath::Vec &acceptor_j, gmath::Vec &acceptor_k) {
      double angles3;
      gmath::Vec tmpA, tmpB;
      tmpA = acceptor_j - donors.pos(i);
      tmpB = acceptor_k - donors.pos(i);
      angles3 = acos((tmpA.dot(tmpB)) / (tmpA.abs() * tmpB.abs()))*180 / M_PI;
      angle_sum += angles3;
      return (angle_sum >= min_angle_sum);
    }

    /**
     * Method that gives back true if the dihedral angle between the atoms i, 
     * the one bound to i, j and k is smaller than the maximal dihedral angle given.
     * @param i
     * @param dihedral
     * @param bound_i
     * @param acceptor_j
     * @param acceptor_k
     * @return 
     */
    bool dihedrals(int i, double &dihedral, gmath::Vec &bound_i, gmath::Vec &acceptor_j, gmath::Vec &acceptor_k) {
      gmath::Vec tmpA, tmpB, tmpC, p1, p2, p3;
      tmpA = bound_i - acceptor_j;
      tmpB = donors.pos(i) - acceptor_k;
      tmpC = acceptor_k - acceptor_j;

      p1 = tmpA.cross(tmpC);
      p2 = tmpB.cross(tmpC);

      dihedral = acos((p1.dot(p2)) / (p1.abs() * p2.abs()))*180 / M_PI;
      p3 = p1.cross(p2);
      if (p3.dot(tmpC) < 0)
        dihedral = -dihedral;
      return (abs(dihedral) <= max_dihedral);
    }
  public:

    /**
     * Constructor, which stores all parameters given from the input file.
     * @param para
     */
    HB3c_calc(HBPara3c para);

    /**
     * Destructor, which closes the timeseries files.
     */
    ~HB3c_calc() {
      timeseriesHB3c.close();
      timeseriesHB3ctot.close();
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
    void calc3();
    /**
     * Method that prints out all 3-centered H-bonds.
     * @param hb2c
     */
    void printstatistics(HB2c_calc &hb2c);
  }; //end class HB3c_calc
}
#endif	/* NEW_3C_H */