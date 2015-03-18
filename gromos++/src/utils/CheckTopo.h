// utils_CheckTopo.h

// Class that runs some basic checks on a molecule topology
#ifndef INCLUDED_CHECKTOPO
#define INCLUDED_CHECKTOPO

#include <vector>
#include <string>

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
   * a class that can un some basic tests on a MoleculeTopology
   *
   * Checks on the bonds, angles, improper dihedrals, exclusions and 
   * charge-groups are carried out
   * @class CheckTopo
   * @author B.C. Oostenbrink
   * @ingroup utils
   */
  class CheckTopo
    {
      /**
       * A pointer to the MoleculeTopology that will be checked.
       */
      const gcore::MoleculeTopology *d_mt;
      /**
       * A vector of strings, which will contain the errors and warnings
       */
      std::vector<std::string> d_error;
      /**
       * The number of digits upto which a charge group should be 
       * an integer
       */
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
       *
       * Function to checks<br>
       * <li> that not more than one bond is defined between two atoms
       * <li> that an atom appears only once in a bond
       */
      int checkBonds();
      /**
       * Check angles
       *
       * Function that checks<br>
       * <li> whether an angle is defined for every atom j which is bonded to
       *      two atoms i and k</li>
       * <li> that not more than one angle is defined for any set of three 
       *      atoms</li>
       * <li> that the atoms defining an angle are bound to the central atom
       * <li> that an atom appears only once in a dihedral
       */
      int checkAngles();
      /**
       * Check Impropers
       *
       * Function that checks<br>
       * <li> whether an improper is defined for every atom i which is bonded
       *      to three other atoms</li>
       * <li> that not more than one improper is defined for any set of four
       *      atoms</li>
       * <li> that all atoms in an improper are bound to each other
       * <li> that an atom appears only once in a dihedral
       */
      int checkImpropers();
      /**
       * Check Dihedrals
       * Function that checks<br>
       * <li> that every atom in a dihedral is bound to the next
       * <li> that an atom appears only once in a dihedral
       */
      int checkDihedrals();
      /**
       * Check CrossDihedrals
       * Function that checks<br>
       * <li> that every atom in a dihedral is bound to the next
       * <li> that an atom appears only once in a dihedral
       */
      int checkCrossDihedrals();
      /**
       * Check exclusions
       *
       * Function that checks<br>
       * <li> That 1,2 and 1,3 neighbours are excluded</li>
       * <li> That excluded atoms are 1,2 or 1,3 neighbours. If 1,4
       *      neighbours are excluded, a warning is produced that this should
       *      be an aromatic system.</li>
       * <li> That 1,4 neighbours are 1,4-excluded (not done if they were 
       *      already excluded</li>
       * <li> That 1,4 excluded atoms are 1,4 neighbours</li>
       */
      int checkExclusions();
      /**
       * a function that sets the precision up to how many digits the charge
       * groups are expected to be an integer.
       */
      void setChargePrecision(int i);
      /**
       * Check chargeGroups
       * 
       * Function that checks whether all charge groups have an integer charge
       */
      int checkChargeGroups();
      /**
       * Function that checks whether the last atom:<br>
       * <li> has no exclusions
       * <li> has no 1,4 exclusions
       * <li> is the end of a charge group
       */
      int checkLastAtom();
      /**
       * Function that calls all checks consecutively.
       */
      int checkAll();
      /**
       * Clear all errors.
       * All errors that were collected this far are discarded.
       */
      void clearErrors();
      
      /**
       * number of errors that were collected
       */
      int numErrors();
      /**
       * error message number i
       */
      std::string error(int i);
      /**
       * the charge precision
       */
      int chargePrecision();
    };
}

#endif
