// args_BoundaryParser.h

#ifndef INCLUDED_ARGS_BOUNDARYPARSER
#define INCLUDED_ARGS_BOUNDARYPARSER

#include <string>

namespace gcore{
  class System;
}

namespace bound{
  class Boundary;
}

namespace args{
  class Arguments;

/**
 * Class BoundaryParser
 * Purpose: Parse boundary type from args.
 * 
 *
 * Description:
 * This class is used to parse boundary types from args. 
 * By default it will use args["pbc"].
 *
 * 
 * @class BoundaryParser
 * @version $Date: Mon Jul 15 14:17:41 MEST 2002
 * @author  R. Buergi
 * @author  M.A. Kastenholz
 * @ingroup args
 * @sa args::Arguments
 * @sa args::GatherParser
 */
  class BoundaryParser{
    
    // not implemented
    BoundaryParser(const BoundaryParser&);
    BoundaryParser &operator=(const BoundaryParser&);

  public:
/**
 * BoundaryParser constructor.
 * Details.
 */
    BoundaryParser(){}
/**
 * BoundaryParser destructor.
 * Details.
 */
    ~BoundaryParser(){}
/** 
 * Constructs the class and returns a pointer to a Boundary.
 * the method parses the input from args.
 * @param sys takes a gcore::system as argument
 * @param args Arguments from the input line.
 * @param str name of the argument string (default "pbc")
 * @return bound::Boundary *boundary pointer to a boudary.
 * Details.
 */
    static bound::Boundary *boundary(gcore::System &sys, const Arguments &args, const std::string &str = "pbc");
    static bound::Boundary *boundary(gcore::System &sys);
  };

}

#endif
