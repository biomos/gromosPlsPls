// gmath_Stat.cc

#ifndef INCLUDED_GMATH_VEC
#include "Vec.h"
#endif

namespace gmath
{
  inline double sqrt(double d)
  {
    return ::sqrt(d);
  }
  
  inline gmath::Vec sqrt(gmath::Vec const & v)
  {
    return gmath::Vec(::sqrt(v[0]),::sqrt(v[1]), ::sqrt(v[2]));
  }

  inline gmath::Vec operator/(gmath::Vec const &v1, gmath::Vec const &v2)
  {
    return gmath::Vec(v1[0]/v2[0], v1[1]/v2[1], v1[2]/v2[2]);
  }
  
  template<typename T>
  Stat<T>::Stat()
    : d_counter(0),
      d_ave(T()),
      d_msd(T()),
      d_ee(T()),
      d_avedone(false),
      d_msddone(false),
      d_eedone(false),
      d_distdone(false),
      d_dist(0, 1, 1)
  { 
  }
  
  template<typename T>
  void Stat<T>::addval(T val)
  {
    d_vals.push_back(val);
    d_counter++;
    d_avedone=false;
    d_msddone=false;
    d_eedone=false;
    if(d_distdone) d_dist.add(val);
  }
  
  template<typename T>
  T Stat<T>::ave()const
  {
    if(!d_avedone){
      d_ave = this->subave(0,d_counter);
      d_avedone=1;
    }
    return d_ave;
  }
  
  template<typename T>
  T Stat<T>::subave(int b, int e)const
  {
    T ave = 0;
    
    //calculate the average
    for(int i=b;i<e;i++){
      ave += d_vals[i];
    }
    return ave/(e-b);
  }
  
  template<typename T>
  T Stat<T>::msd()const
  {
    if(!d_msddone){
      T sum=0, ssum=0;
      for(int i=0; i<d_counter; i++){
	sum+=d_vals[i];
	ssum+=d_vals[i]*d_vals[i];
      }
      sum/=d_counter;
      ssum/=d_counter;
      d_msd = ssum - sum*sum;
    }
    return d_msd;
  }

  template<typename T>
  T Stat<T>::rmsd()const
  {
    return sqrt(this->msd());
  }

  template<typename T>
  T Stat<T>::min()const
  {
    T m=d_vals[0];
    for(int i=1; i<d_counter; ++i)
      if(d_vals[i] < m) m = d_vals[i];
    return m;
  }
  
  template<typename T>
  T Stat<T>::max()const
  {
    T m = d_vals[0];
    for(int i=1; i<d_counter; ++i)
      if(d_vals[i] > m) m=d_vals[i];
    return m;
  }
  
  template<typename T>
  T Stat<T>::ee()const
  {
    if(!d_eedone){
      // first prepare the blocks
      double blksz=50;
      int old=2;
      while(4*blksz<d_counter){
	d_blocksize.push_back(int(blksz));
	old=int(blksz);
	while(old==int(blksz)) blksz = blksz*1.07177;
      }
      
      int Nblocks=d_blocksize.size();
      T rmsd2, ave=0;
      T runave=this->ave();
      T runrmsd=this->rmsd();
      std::vector<T> fit(Nblocks), x(Nblocks);
      
      for(int j=0; j<Nblocks; j++){
	int Nblcki=d_counter/d_blocksize[j];
	
	// The rmsd of the property we are interested in, weighted by the
	// average energy of the blocks
	rmsd2=0;
	for(int i=0; i<Nblcki; i++){
	  ave = this->subave(i*d_blocksize[j],(i+1)*d_blocksize[j]);
	  rmsd2+=(ave-runave)*(ave-runave);
	}
	rmsd2/=Nblcki;
	fit[j]=(d_blocksize[j]*rmsd2) / runrmsd / runrmsd;
	x[j]=1.0/d_blocksize[j];
	
      }
      T sx=0, sf=0,sfx=0,sxx=0;
      for(int i=0; i<Nblocks;i++){
	sx+=x[i];
	sf+=fit[i];
	sfx+=x[i]*fit[i];
	sxx+=x[i]*x[i];
      }
      
      T a, b;
      a=(sf*sx/Nblocks-sfx)/(sx*sx/Nblocks-sxx);
      b = (sf - a*sx)/Nblocks;

      d_ee=sqrt(b/d_counter)*runrmsd;
    }
    
    return d_ee;
  }

  template<typename T>
  std::vector<T> const & Stat<T>::data(){
    return d_vals;
  }

  template<typename T>
  gmath::Distribution const & Stat<T>::distribution()const
  {
    if(d_distdone)
      return d_dist;
    else
      throw Distribution::Exception("call dist_init first");
  }
  
  template<typename T>
  gmath::Distribution const & Stat<T>::dist_init(double lower, double upper, int nsteps)
  {
    d_dist = gmath::Distribution(lower, upper, nsteps);
    //put all values in it
    for(int i=0; i<d_counter; i++)
      d_dist.add(d_vals[i]);
    d_distdone=1;
    
    //return the distribution
    return d_dist;
  }

  template<typename T>
  void Stat<T>::subtract_average()
  {
    double ave=this->ave();
    for(int i=0; i<d_counter; i++)
      d_vals[i]-=ave;
    d_avedone=0;
  }
}

