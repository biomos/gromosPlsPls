#ifndef INCLUDED_UTILS_NOEDISTANCE
#define INCLUDED_UTILS_NOEDISTANCE

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gcore{
  class System;
}


namespace utils{
  class NoeDistance_i;
  
  class NoeDistance{
    NoeDistance_i *d_this;

    // not implemented
    NoeDistance();
    NoeDistance(const NoeDistance &);
    NoeDistance &operator=(const NoeDistance&);
    
  
  public:
    NoeDistance(const VirtualAtom &at1, const VirtualAtom &at2);
  
    // Accessors
    //    inline vector<double> &dist(){return dist_;}
    //    inline const double &dist(const unsigned int i)const{return dist_[i];}
    //    inline unsigned int num_dist(){return dist_.size();}
    //    inline Dist_Res &NoeDistance(unsigned int i){return NoeDistance_[i];}
    //    inline unsigned int num_NoeDistance(){return NoeDistance_.size();}
    
    // Distance correction
    //    const double corr_dist(const int i)const;
    // info string
    //    const string info()const;
    
    struct Exception: public gromos::Exception{
      Exception(const string &str): gromos::Exception("NoeDistance", str){}
    };
    
  };
  
  
} /* namespace */

#endif
