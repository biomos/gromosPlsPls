// gmath_Distribution.cc

#include "Distribution.h"
#include <cmath>
#include <vector>
#include <iomanip>

using gmath::Distribution;

using namespace std;

namespace gmath{
  
Distribution::Distribution(double begin, double end, int nsteps):
  d_count(nsteps)
{
  if(begin>=end) 
    throw Distribution::Exception("Upper boundary should be higher than lower");
  if(nsteps<1) 
    throw Distribution::Exception("You need at least one step");
  d_step=(end-begin)/(nsteps);
  
  for(int i=0;i<nsteps;i++){
    d_count[i]=0;
    
  }
  d_nsteps=nsteps;
  d_begin=begin;
  d_end=end;
  d_sum=0.0;
  d_num=0;
}

void Distribution::write(std::ostream &os)const
{
  
  for(int i=0;i<d_nsteps;i++)
    os << setw(8) << d_begin+(i+0.5)*d_step << "\t" 
       << setw(5) << d_count[i] << endl;
}

void Distribution::write_normalized(std::ostream &os)const
{
  
  for(int i=0;i<d_nsteps;i++)
    os << setw(8) << d_begin+(i+0.5)*d_step << "\t" 
       << setw(5) << double(d_count[i]) / nVal() << endl;
}
  
double Distribution::add(const double value)
{
  if(value>=d_begin&&value<=d_end){
     
    int q=int((value-d_begin)/d_step);
   
    this->d_count[q]++;
    this->d_sum+=value;
    this->d_num++;
    return value;
  }
  
  return value+1;
}

int Distribution::getbin(const double value)
{
  if(value>=d_begin&&value<=d_end){
     
    int q=int((value-d_begin)/d_step);
    return q;
  }
  return -1;
}

bool Distribution::inrange(const double value) {
  if(value>=d_begin&&value<=d_end) return true;
  else return false;
}

double Distribution::rmsd()const
{
  double sumdiff=0;
  double avr=this->ave();
  for(int i=0;i<d_nsteps;i++){
    double diff=avr - (d_begin+(i+0.5)*d_step);
    sumdiff+=d_count[i]*diff*diff;
  }
  return sqrt(sumdiff/d_num);
}
 
}

