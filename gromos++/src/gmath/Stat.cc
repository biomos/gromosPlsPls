// gmath_Stat.cc

#include "Stat.h"
#include "Distribution.h"
#include <vector>
#include <cmath>
#include <new>

namespace gmath
{
  
  stat::stat()
  { 
    d_counter=0;
    d_ave=0;
    d_msd=0;
    d_ee=0;
    d_avedone=0;
    d_msddone=0;
    d_eedone=0;
    d_dist= new Distribution(0,1,1);
    d_distdone=0;
  }
  
  void stat::addval(double val)
  {
    d_vals.push_back(val);
    d_counter++;
    d_avedone=0;
    d_msddone=0;
    d_eedone=0;
    if(d_distdone) d_dist->add(val);
  }
  
  double stat::ave()
  {
    if(!d_avedone){
      d_ave = this->subave(0,d_counter);
      d_avedone=1;
    }
    return d_ave;
  }
  
  double stat::subave(int b, int e)
  {
    double ave=0;
    
    //calculate the average
    for(int i=b;i<e;i++){
      ave+=d_vals[i];
    }
    return ave/(e-b);
  }
  

  double stat::msd()
  {
    if(!d_msddone){
      double sum=0, ssum=0;
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

  double stat::rmsd()
    {
      return sqrt(this->msd());
    }

  double stat::min()
  {
    double m=d_vals[0];
    for(int i=1; i<d_counter; i++)
      if(d_vals[i]<m)m=d_vals[i];
    return m;
  }
  
  double stat::max()
  {
    double m=d_vals[0];
    for(int i=1; i<d_counter; i++)
      if(d_vals[i]>m) m=d_vals[i];
    return m;
  }
  
  double stat::ee()
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
      double rmsd2, ave=0;
      double runave=this->ave();
      double runrmsd=this->rmsd();
      double fit[Nblocks], x[Nblocks];
      
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
	fit[j]=d_blocksize[j]*rmsd2/runrmsd/runrmsd;
	x[j]=1.0/d_blocksize[j];
	
      }
      double sx=0, sf=0,sfx=0,sxx=0;
      for(int i=0; i<Nblocks;i++){
	sx+=x[i];
	sf+=fit[i];
	sfx+=x[i]*fit[i];
	sxx+=x[i]*x[i];
      }
      
      double a, b;
      a=(sf*sx/Nblocks-sfx)/(sx*sx/Nblocks-sxx);
      b = (sf - a*sx)/Nblocks;
      //  for(int i=0; i<Nblocks; i++)
      //    cout << x[i] << "\t" << fit[i] << "\t" << a*x[i]+b << endl;
      
      d_ee=sqrt(b/d_counter)*runrmsd;
    }
    
    return d_ee;
  }

  std::vector<double>* stat::data(){
    return &d_vals;
  }

  gmath::Distribution *stat::distribution()
  {
    if(d_distdone)
      return d_dist;
    else
      throw Distribution::Exception("call dist_init first");
  }
  
  gmath::Distribution *stat::dist_init(double lower, double upper, int nsteps)
  {
    *d_dist = gmath::Distribution(lower, upper, nsteps);
    //put all values in it
    for(int i=0; i<d_counter; i++)
      d_dist->add(d_vals[i]);
    d_distdone=1;
    
    //return the distribution
    return d_dist;
  }
  void stat::substract_average()
  {
    double ave=this->ave();
    for(int i=0; i<d_counter; i++)
      d_vals[i]-=ave;
    d_avedone=0;
  }
}

