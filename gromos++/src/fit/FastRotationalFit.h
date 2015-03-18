// fit_FastRotationalFit.h

#ifndef INCLUDED_FIT_FASTROTATIONALFIT
#define INCLUDED_FIT_FASTROTATIONALFIT


#include "../gromos/Exception.h"

namespace fit{
  /**
   * Class FastRotationalFit
   * This class performs a rotational fit on a set of coordinates
   *
   * A least squares fitting of one set of atoms is performed relative to 
   * another. 
   *
   * @class FastRotationalFit
   * @author Markus Christen, Chris Oostenbrink
   * @ingroup fit
   * @sa fit::RotationalFit
   * @sa fit::TranslationalFit
   * @sa utils::Rmsd
   */
  class FastRotationalFit{
    /**
     * atoms to fit to
     */
    std::vector<bool> d_fit_spec;
    /**
     * number of atoms to fit to
     */
    int d_fit_num_atoms;
    /**
     * atoms to take into account for rmsd calculation
     */
    std::vector<bool> d_rmsd_spec;
    /**
     * number of atoms used in rmsd calculation
     */
    int d_rmsd_num_atoms;

    /**
     * use Kabsch rotational fit
     */
    bool d_kabsch_fit;
    
  public:
    /**
     * Constructor: no specifications
     */
    FastRotationalFit() 
      : d_fit_spec(0), d_fit_num_atoms(0), d_rmsd_spec(0), d_rmsd_num_atoms(0), d_kabsch_fit(false)  {};
    /**
     * Constructor
     */
    FastRotationalFit(std::vector<bool> fit_spec, std::vector<bool> rmsd_spec)
      : d_fit_spec(fit_spec),
	d_fit_num_atoms(0),
	d_rmsd_spec(rmsd_spec),
	d_rmsd_num_atoms(0),
	d_kabsch_fit(false) {
      for(size_t i=0; i<d_fit_spec.size(); ++i)
	if(d_fit_spec[i]) ++d_fit_num_atoms;
      for(size_t i=0; i<d_rmsd_spec.size(); ++i)
	if(d_rmsd_spec[i]) ++d_rmsd_num_atoms;
    };
    
    /**
     * perform a rotational fit (using atoms specified for fitting) on (reduced) atom positions
     * @param ref (reduced) reference atom positions
     * @param sys (reduced) atom positions to be fitted
     */
    int fit(std::vector<gmath::Vec> const & ref,
	    std::vector<gmath::Vec> &sys)const;

    /**
     * perform a rotational fit (using atoms specified for fitting) on atoms specified by sys 
     * to reference atoms specified by ref.
     * @param ref_spec specify reference atoms
     * @param sys_spec specify fittee atoms
     * @param sys [in/out] fittee atoms
     */
    int fit(utils::AtomSpecifier & ref_spec,
	    utils::AtomSpecifier & sys_spec,
	    gcore::System & sys)const;

    /**
     * calculate rotation matrix to do a rotational fit
     * using the atoms specified for fitting
     * @param rot [out] rotation matrix
     * @param ref reference coordinates
     * @param sys coordinates to fit
     */
    int fit(gmath::Matrix &rot,
	    std::vector<gmath::Vec> const &ref,
	    std::vector<gmath::Vec> const &sys)const;

    /**
     * calculate rotation matrix to do a rotational fit
     * using the atoms specified for fitting
     * uses Kabsch algorithm
     * @param rot [out] rotation matrix
     * @param ref reference coordinates
     * @param sys coordinates to fit
     */
    int kabsch_fit(gmath::Matrix &rot,
		   std::vector<gmath::Vec> const &ref,
		   std::vector<gmath::Vec> const &sys)const;
    
    /**
     * calculate the rmsd (using atoms specified for rmsd) on (reduced) atom positions
     * using the specified rotation matrix on sys.
     */
    double rmsd(gmath::Matrix const & rot,
		std::vector<gmath::Vec> const &ref,
		std::vector<gmath::Vec> const &sys)const;
    
    void set_kabsch_fit(bool b)
    {
      d_kabsch_fit = b;
    }

    /**
     * FastRotationalFit exception
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string &what): 
	gromos::Exception("RotationalFit",what){}
    };
    
  };

}

#endif    
