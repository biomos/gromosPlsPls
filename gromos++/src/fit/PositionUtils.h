// fit_PositionUtils.h

#ifndef INCLUDED_FIT_PositionUtils
#define INCLUDED_FIT_PositionUtils


namespace gmath{
  class Vec;
  class Matrix;
}

namespace gcore{
  class System;
}

namespace fit{
  class Reference;

  class PositionUtils{
    // not implemented
    PositionUtils (const PositionUtils&);
    PositionUtils &operator=(const PositionUtils &);
  public:
    PositionUtils(){}
    ~PositionUtils(){}
    
    static gmath::Vec com(const gcore::System &);
    static gmath::Vec com(const gcore::System &, const Reference &);
    static gmath::Vec cog(const gcore::System &);
    static gmath::Vec cog(const gcore::System &, const Reference &);
  
    static void translate(gcore::System *, const gmath::Vec &);
    static void rotate(gcore::System *, const gmath::Matrix &);
    
    static void shiftToCom(gcore::System *);
    static void shiftToCom(gcore::System *, const Reference &);
    static void shiftToCog(gcore::System *);
    static void shiftToCog(gcore::System *, const Reference &);
  };

}

#endif
