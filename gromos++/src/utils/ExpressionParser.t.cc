
#include <iostream>
#include <map>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <cassert>

#include "../gio/InTopology.h"
#include "AtomSpecifier.h"
#include "../gcore/System.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gio/InG96.h"
#include "../bound/Boundary.h"
#include "../args/Arguments.h"
#include "../args/BoundaryParser.h"
#include "../args/GatherParser.h"

#include "ExpressionParser.h"

int main(int argc, char *argv[])
{
  std::string usage = argv[0];
  usage += "\n";
  usage += "\t@topo   <topology>\n";
  usage += "\t@pbc    <pbc>\n";
  usage += "\t@traj   <trajectory>\n";
  usage += "\t@type   <expression type>\n";
  usage += "\t@expr   <expression>\n";
  
  char *knowns[] = {"topo", "pbc", "traj", "type", "expr" };
  args::Arguments args(argc, argv, 5, knowns, usage);

  try{
    if (args["type"] == "double"){
      std::cout << "parsing 'double' expression" << std::endl;
    
      utils::ExpressionParser<double> ep;
    
      std::map<std::string, double> var;
      var["bla"] = 0.12345;
    
      std::string expr1 = args["expr"];

      try{
	double d = ep.parse_expression(expr1, var);
      
	std::cout << "result = " << d << std::endl;
      }
      catch(std::runtime_error e){
	std::cerr << "runtime error: " << e.what() << std::endl;
      }
    }
  
    else if (args["type"] == "int"){
      std::cout << "parsing 'int' expression\n";
    
      utils::ExpressionParser<int> ep;
    
      std::map<std::string, int> var;
      var["bla"] = 12345;
    
      std::string expr1 = args["expr"];
    
      try{
	int d = ep.parse_expression(expr1, var);
      
	std::cout << "result = " << d << std::endl;
      }
      catch(std::runtime_error e){
	std::cerr << "runtime error: " << e.what() << std::endl;
      }
    }
    else{

      std::cout << "parsing 'value' expression\n";

      try{
	gio::InTopology it(args["topo"]);
	gcore::System sys(it.system());
	bound::Boundary *pbc = args::BoundaryParser::boundary(sys, args);

	// parse gather method
	bound::Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

	utils::ExpressionParser<utils::Value> ep(&sys, pbc);
    
	std::map<std::string, utils::Value> var;
	var["bla"] = utils::Value(0.12345);
    
	std::string expr1 = args["expr"];
    
	gio::InG96 ic(args["traj"]);

	// ic.select("ALL");
	ic >> sys;
	(*pbc.*gathmethod)();

	ic.close();

	try{
	  std::cerr << "trying to parse " << expr1 << std::endl;
	  utils::Value v = ep.parse_expression(expr1, var);
      
	  std::cout << "result = " << v << std::endl;
	}
	catch(std::runtime_error e){
	  std::cerr << "runtime error: " << e.what() << std::endl;
	}
      }
      catch (gromos::Exception const & e){
	std::cerr << "Exception:\t";
	std::cerr << e.what() << std::endl;
	return 1;
      }
    }
  }
  catch (gromos::Exception const & e){
    std::cerr << "Exception:\t";
    std::cerr << e.what() << std::endl;
    return 1;
  }

  return 0;
}
