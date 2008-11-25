#ifndef INCLUDED_UTILS_NOE
#define INCLUDED_UTILS_NOE

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gmath {
  class Vec;
}

namespace gcore{
  class System;
}

namespace utils {
  class VirtualAtom;
}


namespace utils{
  class Noe_i;
  /**
   * Class Noe
   * The Noe class stores and analyses Noe information
   *
   * It stores the virtual atoms that define the NOE distance and the 
   * distances
   *
   * @class Noe
   * @author R. Buergi and M.A. Kastenholz
   * @ingroup utils
   * @todo finish the documentation
   */  
  class Noe{

    Noe_i *d_this;

    // not implemented
    Noe();
    Noe(const Noe &);
    Noe &operator=(const Noe&);
    
  
  public:
    /**
     * create an NOE from a line of GROMOS distance restraint specification
     * @param sys the system to create the virtual atoms
     * @param line the line containing the distance restraint specification
     * @param dish carbon-hydrogen distance
     * @param disc carbon-carbon distance
     */
    Noe(gcore::System &sys, const std::string &line, double dish, double disc);
  
    /**
     * the distance corresponding to the NOE.
     * Periodic boundary conditions are not taken into account.
     */
    double distance(int i)const;
    /**
     * the distance vector corresponding to the NOE.
     * Periodic boundary conditions are not taken into account.
     */
    gmath::Vec distanceVec(int i)const;
    
    /**
     * get the virtual atoms
     * @param i the atom i occuring
     * @param ii of the distance ii
     */
    const utils::VirtualAtom & getAtom(int i, int ii) const;
   
    /**
     * the reference distance including the correction
     */
    double correctedReference(int i)const;

    /**
     * a string containing the GROMOS distance restraint specification line
     */
    std::string distRes(int i)const;

    /**
     * a string containing information about the NOE including
     * residue number, residue name, atom and molecule number and type.
     */
    std::string info(int i)const;


    /**
     * the number of distances
     */
    int numDistances()const;
    /**
     * the number of references
     */
    int numReferences()const;
    /**
     * the reference length of the NOE.
     * @param i index of the reference length (0)
     */
    double reference(int i)const;
    /**
     * the correction length for type
     * @param type the virtual atom type
     */
    double correction(int type);

    /**
     * set the correction length for a type
     * @param type virtual atom type
     * @param correction correction length
     */
    void setcorrection(int type, double correction);
    
    struct Exception: public gromos::Exception{
      Exception(const std::string &str): gromos::Exception("Noe", str){}
    };
    
  };

    

} /* namespace */

#endif
