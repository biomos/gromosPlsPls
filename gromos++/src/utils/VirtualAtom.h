// utils_VirtualAtom.h

#ifndef INCLUDED_UTILS_VIRTUALATOM
#define INCLUDED_UTILS_VIRTUALATOM

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif
#ifndef INCLUDED_GMATH_VEC
#include "../gmath/Vec.h"
#endif

namespace gcore{
  class System;
}

namespace utils{

  class VirtualAtom_i;
  /**
   * Class VirtualAtom
   * This class contains information about virtual atoms
   *
   * From this information one should be able to calculate the coordinates
   * for any hydrogen that is not actually there
   *
   * @class VirtualAtom
   * @author R. Buergi and M.K. Kastenholz
   * @ingroup utils
   * @todo finish documentation
   */
  class VirtualAtom{
    VirtualAtom_i *d_this;

    // not implemented
    VirtualAtom();
    VirtualAtom &operator=(const VirtualAtom&);
  
  public:
    /* Create a virtual atom out of atom "atom" of molecule "mol" 
       of System "sys". (numbering from 0!): 
       type indicates the type of the virtual atom as in 
       table 2.6.4.1. In case of a stereospecific 
       CH2, indicate the orientation as 0 or 1.*/
    VirtualAtom(const gcore::System &sys, int mol, int atom, 
		int type, int orientation=0);
    VirtualAtom(const VirtualAtom&);
    ~VirtualAtom();

    // calculates the position of the virtual atom
    gmath::Vec pos() const;

    // sets carbon-hydrogen distance (0.1 by default)
    static void setDish(double dish);
    // sets carbon-carbon distance (0.153 by default)
    static void setDisc(double disc);
    

    // Accessors
    // type
    int type()const;
    // molecule number
    int mol()const;
    // orientation for type 4 CH1
    int orientation()const;
    
    // returns configuration of distance restraint
    int operator[](int i)const;

    // Exception
    struct Exception: public gromos::Exception{
      Exception(const string& what) : gromos::Exception("VirtualAtom", what){}
    };
    
    
  };
}
#endif
