// args_BoundaryParser.cc

#include "BoundaryParser.h"
#include "Arguments.h"
#include "../bound/Vacuum.h"
#include "../bound/TruncOct.h"
#include "../bound/RectBox.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include <sstream>

using namespace args;
using namespace bound;

bound::Boundary *BoundaryParser::boundary(gcore::System &sys, const Arguments &args){
  Boundary *pbc;

  try{

    Arguments::const_iterator it=args.lower_bound("pbc");
    if(it == args.upper_bound("pbc"))
      throw Arguments::Exception("");

//    if (it->second.size() > 1)
  //    throw gromos::Exception("Boundary", args["pbc"] + 
//			      " unknown. Known boundaries are r, t and v");
    
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
      throw gromos::Exception("Boundary", args["pbc"] + 
			      " unknown. Known boundaries are r, t and v");
    }

    ++it;
    /*    for(int i=0;
	it!=args.upper_bound("pbc"); ++it, ++i){
      int at = atoi(it->second.c_str());
      assert(at>0);
      at--;
      int mol=0, atNum=0;
      
      // parse into mol and atom rather than high atom nr.
      while(at > (atNum+=sys.mol(mol).numAtoms())) ++mol;
      at-=atNum-sys.mol(mol).numAtoms();
      
      
      pbc->setReference(i,sys.mol(mol).pos(at));

      } */



  }
  catch(Arguments::Exception &e){
    pbc = new Vacuum(&sys);
  }


  return pbc;
}
