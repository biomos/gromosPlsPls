// args_GatherParser.h

#ifndef INCLUDED_ARGS_GATHERPARSER
#define INCLUDED_ARGS_GATHERPARSER

#include "../bound/Boundary.h"


namespace gcore{
  class System;
}

namespace bound{
  class Boundary;
}

namespace args{
  class Arguments;


/**
 * Class GatherParser
 * Purpose: Parse gathering methods from args["pbc"].
 * 
 *
 * Description:
 * This class is used to parse gathering methods from args["pbc"].
 * By default (when no arguments are given) it will use the gather-gromos 
 * gathering method.
 *
 * 
 * @class GatherParser
 * @version $Date: Mon Jul 15 14:17:41 MEST 2002
 * @author  M.A. Kastenholz
 * @sa args::Arguments
 * @sa args::BoundaryParser
 */
  class GatherParser{
    
    // not implemented
    GatherParser(const GatherParser&);
    GatherParser &operator=(const GatherParser&);

  public:
/**
 * GatherParser constructor.
 * Details.
 */
    GatherParser(){}
/**
 * GatherParser destructor.
 * Details.
 */
    ~GatherParser(){}
/** 
 * Constructs the class and returns a member pointer to a Boundary.
 * Method parse parses the input from args["pbc"].
 * @param args Arguments from the input line.
 * @return bound::Boundary::MemPtr Member pointer to the Boundary class.
 * Details.
 */
  static bound::Boundary::MemPtr parse(const Arguments &args);
  };

}






#endif
