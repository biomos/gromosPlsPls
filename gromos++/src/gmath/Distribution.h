// gmath_Distribution

#ifndef INCLUDED_GMATH_DISTRIBUTION
#define INCLUDED_GMATH_DISTRIBUTION

#include <iostream>
#include <vector>
#include "../gromos/Exception.h"

namespace gmath{

class Distribution{
  double d_begin, d_end, d_step, d_sum;
  int d_nsteps, d_num;
  std::vector<int> d_count;

 public:
  Distribution(double begin=0, double end=1, int nsteps=100);
  ~Distribution(){}
  // Methods
  double add(const double value);
  void write(std::ostream &os)const;

  double ave()const;
  
  double rmsd()const;

  // Accessors
  int operator[](int i)const;
  
  double value(int i)const;
  int nVal()const;
  int nSteps()const;
  
  // Exceptions
  struct Exception: public gromos::Exception{
    Exception(const string& what): gromos::Exception("Distribution", what){}
  };

};
// inline functions & free operators
  inline double Distribution::ave()const{
    return d_sum/d_num;
  }
  inline int Distribution::operator[](int i)const{
    return d_count[i];
  }
  inline double Distribution::value(int i)const{
    return d_begin+(i+0.5)*d_step;
  }
  inline int Distribution::nVal()const
    {
      return d_num;
    }
  
  inline int Distribution::nSteps()const
    {
      return d_nsteps;
    }
  
}

#endif







