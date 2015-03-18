// fit_AtomDistances.h

#ifndef INCLUDED_FIT_ATOMDISTANCES
#define INCLUDED_FIT_ATOMDISTANCES

#include "../gromos/Exception.h"

namespace fit{
  /**
   * Class AtomDistances
   * This class calculates the structural difference between to structures
   * based on the interatomic distances
   *
   * @class AtomDistances
   * @author Chris Oostenbrink
   * @ingroup fit
   * @sa fit::FastRotationalFit
   * @sa fit::TranslationalFit
   * @sa utils::Rmsd
   */
  class AtomDistances{
    /**
     * atoms to calculate the distances for
     */
    std::vector<bool> d_dist_spec;
    /**
     * number of atoms used in dist calculation
     */
    int d_dist_num_atoms;

  public:
    /**
     * Constructor: no specifications
     */
    AtomDistances() 
      : d_dist_spec(0), d_dist_num_atoms(0) {};
    /**
     * Constructor
     */
    AtomDistances(std::vector<bool> dist_spec)
      : d_dist_spec(dist_spec),
	d_dist_num_atoms(0){
      for(size_t i=0; i<d_dist_spec.size(); ++i)
	if(d_dist_spec[i]) ++d_dist_num_atoms;
    };
    
    /**
     * calculate the distance (using atoms specified for dist) on (reduced) atom positions
     * using the specified rotation matrix on sys.
     */
    double dist(std::vector<gmath::Vec> const &ref,
		std::vector<gmath::Vec> const &sys)const;
    

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
