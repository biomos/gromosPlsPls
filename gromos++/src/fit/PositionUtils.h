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
  /**
   * Class PositionUtils
   * Class that contains several operations to perform on your system
   *
   * Different functions that change the coordinates of your system are 
   * defined. It only works for the solute!
   *
   * @class PositionUtils
   * @author R. Buergi
   * @ingroup fit
   * @sa Reference
   */
  class PositionUtils{
    // not implemented
    PositionUtils (const PositionUtils&);
    PositionUtils &operator=(const PositionUtils &);
  public:
    /**
     * PositionUtils constructor
     */
    PositionUtils(){}
    /**
     * PositionUtils deconstructor
     */
    ~PositionUtils(){}
    
    /**
     * Method to calculat the centre of mass of your system
     * @return returns a Vec with the centre of mass position
     */
    static gmath::Vec com(const gcore::System &);
    /**
     * Method to calculate the centre of mass of your system, where a 
     * weight is added from the Reference (might be zero for many atoms)
     * @return a vector with the centre of mass position
     */
    static gmath::Vec com(const gcore::System &, const Reference &);
    /**
     * Method to calculate the centre of geometry of your system
     * @return a Vec with the centre of geometry position
     */
    static gmath::Vec cog(const gcore::System &);
    /**
     * Method to calculate the centre of geometry of your system, where
     * a weight is added from the Reference (might be zero for many atoms)
     * @return a Vec with the centre of geometry position
     */
    static gmath::Vec cog(const gcore::System &, const Reference &);
  
    /**
     * Method to move all solute atoms of your system by a Vec
     */
    static void translate(gcore::System *, const gmath::Vec &);
    /**
     * Method to rotate all solute atoms of your system according to a 
     * rotation Matrix
     */
    static void rotate(gcore::System *, const gmath::Matrix &);
    /**
     * Method to calculate the matrix that rotates around 
     * the specified vector with the specified angle
     */
    static gmath::Matrix rotateAround(gmath::Vec v, double a);
    
    /**
     * Method to translate the System in such a way that its centre of
     * mass is at the origin
     */
    static void shiftToCom(gcore::System *);
    /** 
     * Method to translate the System in such a way that its centre of 
     * mass (calculated with the weights from the Reference) is at the 
     * origin
     */
    static void shiftToCom(gcore::System *, const Reference &);
    /**
     * Method to translate the System in such a way that its centre of
     * geometry is at the origin
     */
    static void shiftToCog(gcore::System *);
    /**
     * Method to translate the System in such a way that its centre of 
     * geometry (calculated with the weights from the Reference) is at 
     * the origin
     */
    static void shiftToCog(gcore::System *, const Reference &);
  };

}

#endif
