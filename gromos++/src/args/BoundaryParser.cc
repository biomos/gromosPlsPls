// args_BoundaryParser.cc

#include "BoundaryParser.h"
#include "Arguments.h"
#include "../bound/Vacuum.h"
#include "../bound/TruncOct.h"
#include "../bound/RectBox.h"
#include "../bound/Triclinic.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Box.h"

using namespace args;
using namespace bound;

bound::Boundary *BoundaryParser::boundary(gcore::System &sys,
        const Arguments &args,
        const std::string &str) {
  Boundary *pbc;

  try {

    // get the boundary from a read box, only if no @pbc flag is set in the args
    if (sys.hasBox && sys.box().boxformat() != gcore::Box::box96 && args.count("pbc") < 0) {
      switch (sys.box().ntb()) {
        case gcore::Box::vacuum :
                  pbc = new Vacuum(&sys);
          break;
        case gcore::Box::rectangular :
                  pbc = new RectBox(&sys);
          break;
        case gcore::Box::triclinic :
                  pbc = new Triclinic(&sys);
          break;
        case gcore::Box::truncoct :
                  pbc = new TruncOct(&sys);
          break;
        default:
          throw gromos::Exception("Boundary", "Could not parse the boundary from"
                  " the box.");
      }
    } else { // parse the boundary condition with help of the command line argument @pbc
      Arguments::const_iterator it = args.lower_bound(str);
      if (it == args.upper_bound(str))
        throw Arguments::Exception("");

      switch (it->second[0]) {
        case 'r':
          pbc = new RectBox(&sys);
          break;
        case 't':
          pbc = new TruncOct(&sys);
          break;
        case 'c':
          pbc = new Triclinic(&sys);
          break;
        case 'v':
          pbc = new Vacuum(&sys);
          break;
        default:
          throw gromos::Exception("Boundary", args[str] +
                  " unknown. Known boundaries are r, t and v");
      }

      ++it;
    }
  }  catch (Arguments::Exception &e) {
    pbc = new Vacuum(&sys);
  }

  return pbc;
}




