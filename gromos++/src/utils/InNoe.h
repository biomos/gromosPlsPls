#ifndef INCLUDED_UTILS_INNOE
#define INCLUDED_UTILS_INNOE

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gcore{
  class System;
}


namespace utils{
  class InNoe_i;
  
  class InNoe{
    InNoe_i *d_this;

    // not implemented
    InNoe();
    InNoe(const InNoe &);
    InNoe &operator=(const InNoe&);
    
  
  public:
    /* 
       Constructs an InNoe from a System and a line
       the line is in the format
       atom1 atom2 distance1 [, distance2, ...]
       where atom1,2 can be 
          ##E for an explicit atom
          ##S for a stereospecific CH1 or CH2
          ##N for a non-stereospecific CH2 or CH3 or two rotating CH3 groups
          ##S## for two stereospecific methyl groups (Val, Leu)
	  ##A for an CD or CE in Phe or Tyr
	  default is ##E
        distances can be indicated as many as there can be due to ambiguities
    */
    
    InNoe(const gcore::System &sys, const string &line);
  
    // Accessors
    //    inline vector<double> &dist(){return dist_;}
    //    inline const double &dist(const unsigned int i)const{return dist_[i];}
    //    inline unsigned int num_dist(){return dist_.size();}
    //    inline Dist_Res &InNoe(unsigned int i){return InNoe_[i];}
    //    inline unsigned int num_InNoe(){return InNoe_.size();}
    
    // Distance correction
    //    const double corr_dist(const int i)const;
    // info string
    //    const string info()const;
    
    struct Exception: public gromos::Exception{
      Exception(const string &str): gromos::Exception("InNoe", str){}
    };
    
  };
  
  
} /* namespace */

#endif
