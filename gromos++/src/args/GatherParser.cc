// args_GatherParser.cc

#include "GatherParser.h"
#include "Arguments.h"
#include "../bound/Boundary.h"
#include <sstream>
#include <iostream>

using namespace args;
using namespace bound;

 bound::Boundary::MemPtr GatherParser::parse(const Arguments &args){
   Boundary::MemPtr gathmethod; 

  try{

    Arguments::const_iterator it=args.lower_bound("pbc");
    if(it == args.upper_bound("pbc"))
      throw Arguments::Exception("");

    if (it->second.size() <= 1) {gathmethod = &Boundary::coggather;}
  
    if (it->second == "r" || "t" || "v") { ++it;}   
  
    string gather =  it->second;

    if (gather == "g"){
     gathmethod = &Boundary::gather;
    }
   else if (gather == "ggr"){
   gathmethod = &Boundary::gathergr;
   }
   else if (gather == "cog"){
   gathmethod = &Boundary::coggather;
   }
    else {
   throw gromos::Exception("Gather", args["gather"] + 
			      " unknown. Known boundaries are g, gr and cog");
    }
    
  }
  catch(Arguments::Exception &e){
    gathmethod = &Boundary::gather;
  }


  return gathmethod;
}
