// gmath_STAT

#ifndef INCLUDED_GMATH_STAT
#define INCLUDED_GMATH_STAT

#ifndef INCLUDED_GMATH_DISTRIBUTION
#include "Distribution.h"
#define INCLUDED_GMATH_DISTRIBUTION
#endif

#include <vector>

namespace gmath
{
  /**
   * Class stat
   * A class to perform some basic statistics on a series of numbers
   *
   * This class allows one to store a series of numbers and calculate 
   * the average, rmsd and an error estimate
   *
   * @class stat
   * @author B.C. Oostenbrink
   * @ingroup gmath
   * @todo Include the option to immediately write a distribution for the 
   *       values that are stored
   */
  class stat
    {
      std::vector<int> d_blocksize;
      std::vector<double> d_ener, d_vals;
      int d_counter;
      double d_ave,d_msd, d_ee;
      int d_avedone, d_msddone, d_eedone, d_distdone;
      gmath::Distribution *d_dist;
      
    public:
      /**
       * stat constructor
       */
      stat();
      /**
       * stat deconstructor
       */
      ~stat(){}
      /**
       * Method to add another value to the series
       * @param val the value to add
       */
      void addval(double val);
      /**
       * Method to calculate (or return) the mean square deviation 
       * of the series. Within the class, we keep track of whether 
       * anything has changed since the previous calculation to determine
       * if a new calculation is needed
       * @return The mean square deviation
       */
      double msd();
      /**
       * Method to calculate the square root of the mean square deviation
       * @return The rmsd
       */
      double rmsd();
      /**
       * Method to calculate the average over the series. Internally, we
       * determine whether a new calculation is required or if we can just 
       * return the previously calculated value.
       * @return The average
       */
      double ave();
      /**
       * Method to calculate the average over only part of the series.
       * @param b first index of the series
       * @param e last index of the series. The average is calculated 
       *          for(i=b; i<e; i++)
       * @return The average of this range
       */
      double subave(int b, int e);
      /**
       * Method to calculate an error estimate for the series.
       *
       * The error estimation is based on a method described by Alan and 
       * Tildesley. The series is devided into blocks, for which the average
       * is calculated. The rmsd of these averages is then calculated for 
       * different block sizes. An extrapolation to infinite block size then 
       * gives the error estimate.
       * @return The error estimate
       */
      double ee();
      /**
       * Accessor that returns the stored values
       * @param i the i-th value is returned
       * @return the value that was stored
       */
      double val(int i);
      /**
       * Accessor to return the number of elements that have been stored 
       * so far
       * @return the number of values that are stored in the class
       */
      int n();
      /**
       * Accessor that returns the minimum value of the values stored so 
       * far
       */
      double min();
      /**
       * Accessor that returns the maximum value of the values stored so
       * far
       */
      double max();
      /**
       * Accessor that returns a pointer to the vector containing the data
       */
      std::vector<double> *data();
      /**
       * Accessor that returns a pointer to a Distribution of the data
       *
       * The stat class also contains a Distribution of the data. 
       * Call dist_init to define the bounds of the Distribution.
       * @return a pointer to a gmath::Distribution
       */
      gmath::Distribution *distribution();
      /**
       * Initializes a Distribution
       *
       * The function dist_init rescales the Distribution to the specified
       * bounds and puts all values so far in the distribution.
       * Values that are added later, are also added to the distribution.
       * @param lower lower bound of the distribution
       * @param upper upper bound of the distribution
       * @param nsteps number of grid points for the distribution
       */
      gmath::Distribution *dist_init(double lower, double upper, int nsteps);
      /**
       * Initializes a Distribution
       *
       * The function dist_init rescales the Distribution to use the current
       * minima and maxima as bounds.
       * @param nsteps number of grid points for the distribution
       */
      gmath::Distribution *dist_init(int nsteps);
      /**
       * Substract the average from every data point
       */
      void substract_average();
      
    };

  //inline functions
  inline double stat::val(int i)
    {
      return d_vals[i];
    }
  
  inline int stat::n()
    {
      return d_counter;
    }

  inline gmath::Distribution *stat::dist_init(int nsteps)
    {
      return this->dist_init(this->min(), this->max(), nsteps);
    }
  
};
#endif
