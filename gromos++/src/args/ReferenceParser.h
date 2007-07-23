// args_ReferenceParser.h

#ifndef INCLUDED_ARGS_REFERENCEPARSER
#define INCLUDED_ARGS_REFERENCEPARSER

#include <vector>
#include "../fit/Reference.h"

namespace gcore{
  class System;
}

namespace args{
  class Arguments;
}

namespace fit{
  class Reference;
}

namespace args{
  /**
   * Class ReferenceParser
   * Parses the different inputs for specifying the reference from the 
   * arguments 
   * 
   * The reference parser looks at a set of input flags that can be used
   * to specify the reference atoms (class, atoms). These are required 
   * for a (rotational) fit. The atoms are added to the reference
   * immediately.
   *
   * @class ReferenceParser
   * @author V. Kraeutler
   * @ingroup args
   * @sa args::Arguments
   */
  class ReferenceParser{
    
    std::vector<int> myMols;
    gcore::System &mySys;
    const args::Arguments &myArgs;
    fit::Reference &myRef;
    bool myAdded;

  public:

    /**
     * ReferenceParser constructor
     *
     * @param sys The system on which the reference is based
     * @param args All the arguments, the parser only looks at the 
     *             the interesting ones
     * @param ref The reference to which the atoms need to be added
     */
    ReferenceParser(
      gcore::System &sys,
      const args::Arguments &args,
      fit::Reference &ref);
    /**
     * Function that is only needed internally
     */
    void add_ref();
    /**
     * Function that is only needed internally
     */
    void getMolecules();
    /**
     * Function that is only needed internally
     */
    void classes();
    /**
     * Function that is only needed internally
     */
    void atoms();
  };

}
#endif
