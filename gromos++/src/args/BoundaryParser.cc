// args_BoundaryParser.cc

#include "BoundaryParser.h"
#include "Arguments.h"
#include "../bound/Vacuum.h"
#include "../bound/TruncOct.h"
#include "../bound/RectBox.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"

using namespace args;
using namespace bound;

bound::Boundary *BoundaryParser::boundary(gcore::System &sys, 
					  const Arguments &args, 
					  const std::string &str){
  Boundary *pbc;

  try{

    Arguments::const_iterator it=args.lower_bound(str);
    if(it == args.upper_bound(str))
      throw Arguments::Exception("");

    switch(it->second[0]){
    case 'r':
      pbc=new RectBox(&sys);
      break;
    case 't':
      pbc=new TruncOct(&sys);
      break;
    case 'v':
      pbc=new Vacuum(&sys);
      break;
    default:
      throw gromos::Exception("Boundary", args[str] + 
			      " unknown. Known boundaries are r, t and v");
    }

    ++it;
  }
  catch(Arguments::Exception &e){
    pbc = new Vacuum(&sys);
  }

  return pbc;
}




