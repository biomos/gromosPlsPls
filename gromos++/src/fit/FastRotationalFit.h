// fit_FastRotationalFit.h

#ifndef INCLUDED_FIT_FASTROTATIONALFIT
#define INCLUDED_FIT_FASTROTATIONALFIT


#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace fit{
  /**
   * Class FastRotationalFit
   * This class performs a rotational fit on a set of coordinates
   *
   * A least squares fitting of one set of atoms is performed relative to 
   * another. 
   *
   * @class RotationalFit
   * @author Markus Christen, Chris Oostenbrink
   * @ingroup fit
   * @sa fit::RotationalFit
   * @sa fit::TranslationalFit
   * @sa utils::Rmsd
   */
  class FastRotationalFit{
    std::vector<bool> d_fit_spec;
    int d_fit_num_atoms;
    std::vector<bool> d_rmsd_spec;
    int d_rmsd_num_atoms;
    
  public:
    FastRotationalFit() :d_fit_spec(0), d_fit_num_atoms(0), d_rmsd_spec(0), d_rmsd_num_atoms(0)  {};
    FastRotationalFit(std::vector<bool> fit_spec, std::vector<bool> rmsd_spec):
      d_fit_spec(fit_spec),
      d_fit_num_atoms(0),
      d_rmsd_spec(rmsd_spec),
      d_rmsd_num_atoms(0) {
    for(size_t i=0; i<d_fit_spec.size(); ++i)
      if(d_fit_spec[i]) ++d_fit_num_atoms;
    for(size_t i=0; i<d_rmsd_spec.size(); ++i)
      if(d_rmsd_spec[i]) ++d_rmsd_num_atoms;
    };
    
    
    
    
    int fit(gmath::Matrix &rot, std::vector<gmath::Vec> const &ref, std::vector<gmath::Vec> const &sys)const;
    int fit(std::vector<gmath::Vec> const & ref, std::vector<gmath::Vec> &sys)const;
    
    double rmsd(gmath::Matrix const & rot, std::vector<gmath::Vec> const &ref, std::vector<gmath::Vec> const &sys)const;
    
    struct Exception: public gromos::Exception{
      Exception(const std::string &what): 
	gromos::Exception("RotationalFit",what){}
    };
    
  };

}

#endif    
