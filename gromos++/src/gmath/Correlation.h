// gmath_correlation

#ifndef INCLUDED_GMATH_CORRELATION
#define INCLUDED_GMATH_CORRELATION

#include "Vec.h"
#include "Stat.h"
#include <vector>

namespace gmath
{
  /**
   * Class correlation
   * A class to calculate time correlation functions
   *
   * This class can calculate almost any kind of correlation function
   * between two timeseries of scalars or vectors. Data should be provided
   * by either two (or one) vectors of double, statistic classes or vectors
   * gmathh:Vec for vector correlation functions
   *
   * @class Correlation
   * @author B.C. Oostenbrink
   * @ingroup gmath
   */
  class Correlation
    {
      std::vector<double> d_f;
      const std::vector<double> *d_a;
      const std::vector<double> *d_b;
      std::vector<gmath::Vec> *d_va;
      std::vector<gmath::Vec> *d_vb;
      bool d_vec;
      bool d_calc;
    public:
      /**
       * correlation constructor for two vectors of double
       */
      Correlation(std::vector<double> &a, std::vector<double> &b);
      /**
       * correlation constructor for two vectors of Vec
       */
      Correlation(std::vector<gmath::Vec> &a, std::vector<gmath::Vec> &b);
      /**
       * correlation constructor for two stat-classes
       */
      Correlation(gmath::Stat<double>& a, gmath::Stat<double>& b);
      /**
       * correlation deconstructor
       */
      ~Correlation(){}
       /**
	* method to calculate correlation function with the use of fft
	* 
	* is only applicable to the time correlation function of two scalars
	* which is defined by the product of the two:
	* C(t) = <A(T)B(T+t)>
	*/
      void calc_fft();
      /**
       * method to calculate correlation function with calculating the direct
       * product (slower than with fft)
       *
       * is applicable to the time correlation function of both scalars and
       * vectors. Uses the (dot)product of the two entities.
       * C(t) = <A(T).B(T+t)>
       */
      void calc_direct();
      /**
       * method to calculate the correlation function according to a
       * specified expression. This uses a direct method and will be
       * versatile but slow
       *
       * C(t) = <f(A(T), B(T+t))> 
       * @param string s  An expression according to gmath::expression
       */
      void calc_expression(std::string s);
      /**
       * Accessor to the i-th element of the time correlation function
       */
      double operator[](int i);
      /**
       * Accessor to the size of the correlation function
       */
      int size();
      /**
       * Method to calculate the power spectrum of the correlation function
       * using fft
       * @param vector<double> w will be returned with the frequencies
       * @param vector<double> s will be returned with the intensity at these 
       *                         frequencies
       * @param double dt time step for data (and correlation function)
       * @param double frac determines which fraction of the correlation 
       *                    function will be taken into account 
       *                    (noise reduction)
       */
      void spectrum(std::vector<double> &w, std::vector<double> &s, double dt,
		    double frac=1.0);


    };
}
#endif
