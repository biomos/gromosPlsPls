/**
 * @file TranslationalFit.h
 */

#ifndef INCLUDED_FIT_TRANSLATIONALFIT
#define INCLUDED_FIT_TRANSLATIONALFIT


namespace gcore{
  class System;
}

namespace gmath{
  class Vec;
}

namespace fit{
  class TranslationalFit_i;
  class Reference;

  enum centre_enum { cog, com };

  /**
   * Class TranslationalFit
   * A class that performs a translational fit of one system with 
   * respect to a reference
   *
   * @class TranslationalFit
   * @author R. Buergi
   * @ingroup fit
   * @sa RotationalFit
   * @sa Reference
   */
  class TranslationalFit{
    TranslationalFit_i *d_this;

    // not implemented
    TranslationalFit();
    TranslationalFit(const TranslationalFit &);
    TranslationalFit &operator=(const TranslationalFit &);

  public:
    // construct Reference for reference molecule, then
    // construct the TranslationalFit.
    /**
     * TranslationalFit constructor.
     * @arg Reference reference to which your system will be fitted
     * @arg centre cog or com : centre of geometry or centre of mass.
     *      default is centre of geometry
     * 
     */
    TranslationalFit(Reference *ref, centre_enum centre = fit::cog);
    /**
     * TranslationalFit deconstructor
     */
    ~TranslationalFit();
    
    //    void fitToCom(gcore::System *)const;
    /**
     * Method that fits the System to the Reference
     */
    void fit(gcore::System *)const;
    /**
     * Accessor that returns the centre of mass of the system including the 
     * weight defined in the Reference
     */
    const gmath::Vec &com()const;
    /**
     * Accessor that returns the centre of geometry of the system, including
     * the weights defined in the Reference
     */
    const gmath::Vec &cog()const;
    
  };

}

#endif    
