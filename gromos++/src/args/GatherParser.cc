// args_GatherParser.cc

#include "GatherParser.h"
#include "Arguments.h"
#include "../bound/Boundary.h"

using namespace args;
using namespace bound;

bound::Boundary::MemPtr GatherParser::parse(const Arguments &args,const std::string &str){
  Boundary::MemPtr gathmethod; 

  try{

    Arguments::const_iterator it=args.lower_bound(str);
    if(it == args.upper_bound(str))
      throw Arguments::Exception("");
    ++it;  

    if(it == args.upper_bound(str)) {
      gathmethod = &Boundary::coggather;}
 
    else {
      std::string gather =  it->second;
   
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
	throw gromos::Exception("Gather", args[str] + 
				" unknown. Known gathering methods are g, ggr and cog");
      }
    }  
  }
  catch(Arguments::Exception &e){
    gathmethod = &Boundary::coggather; 
  }


  return gathmethod;
}
