
#include <iostream>
#include <map>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <cassert>

#include "../gcore/System.h"
#include "../bound/Boundary.h"

#include "ExpressionParser.h"

int main(int argc, char *argv[])
{
  utils::ExpressionParser<double> ep;

  std::map<std::string, double> var;
  var["bla"] = 0.12345;
  
  if (argc < 2){
    std::cout << "usage: " << argv[0] << " expr\n\n";
    return 1;
  }
    
  std::string expr1 = argv[1];

  try{
    double d = ep.parse_expression(expr1, var);
  
    std::cout << "result = " << d << std::endl;
  }
  catch(std::runtime_error e){
    std::cerr << "runtime error: " << e.what() << std::endl;
  }
  
  return 0;
}
