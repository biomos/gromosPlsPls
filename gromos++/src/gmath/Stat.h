// gmath_STAT

#ifndef INCLUDED_GMATH_STAT
#define INCLUDED_GMATH_STAT

#include <vector>

namespace gmath
{
  
  class stat
    {
      std::vector<int> d_blocksize;
      std::vector<double> d_ener, d_vals;
      int d_counter;
      double d_ave,d_msd, d_ee;
      int d_avedone, d_msddone, d_eedone;
      
    public:
      stat();
      ~stat(){}
      void addval(double val);
      double msd();
      double rmsd();
      double ave();
      double subave(int b, int e);
      double ee();
      double val(int i);
      int n();

      
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

  
};
#endif
