// utils_CheckTopo.h

// Class that runs some basic checks on a molecule topology
#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif
#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore
{
  class System;
  class MoleculeTopology;
  class AtomTopology;
}
namespace utils
{
  /**
   * Class CheckTopo
   * @todo finish documentation
   */
  class CheckTopo
    {
      const gcore::MoleculeTopology *d_mt;
      std::vector<std::string> d_error;
      int d_chargePrecision;
    public:
      /**
       * Constructor
       */
      CheckTopo(const gcore::MoleculeTopology &mt) { 
	d_mt = &mt; 
	d_chargePrecision = 5;
      };
      /**
       * Constructor from a system and a molecule number
       */
      CheckTopo(const gcore::System &sys, int m);
      
      /**
       * deconstructor
       */
      ~CheckTopo(){};
	
      /**
       * Check bonds
       */
      int checkBonds();
      /**
       * Check angles
       */
      int checkAngles();
      /**
       * Check Impropers
       */
      int checkImpropers();
      /**
       * Check exclusions
       */
      int checkExclusions();
      /**
       * set the charge precision
       */
      void setChargePrecision(int i);
      /**
       * Check chargeGroups
       */
      int checkChargeGroups();
      /**
       * Check everything
       */
      int checkAll();
      /**
       * Clear all errors
       */
      void clearErrors();
      
      /**
       * number of errors
       */
      int numErrors();
      /**
       * error message
       */
      std::string error(int i);
      /**
       * the charge precision
       */
      int chargePrecision();
    };
}


      
  
