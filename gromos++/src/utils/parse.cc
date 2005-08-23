
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <sstream>
#include <iostream>

#include "parse.h"
#include "ExpressionParser.h"

namespace utils
{
  
  std::string::size_type find_par
  (
   std::string s,
   char c,
   std::string::size_type it
   )
  {
    // recognized brackets: (, [, {, <
    int level = 0;
    
    for( ; it < s.length(); ++it){
      
      switch(s[it]){
	case '(' :
	case '[' :
	case '{' :
	case '<' :
	  ++level;
	  break;
	case ')' :
	case '}' :
	case ']' :
	case '>' :
	  --level;
	  break;
	  
	default :
	  ;
      }
      
      if (s[it] == c && level == 0) return it;
      
    }
    
    return std::string::npos;
  }
  
  std::string::size_type find_matching_bracket
  (
   std::string s,
   char bra,
   std::string::size_type it
   )
  {
    char ket;
    if (bra == '(') ket = ')';
    else if (bra == '[') ket = ']';
    else if (bra == '{') ket = '}';
    else if (bra == '<') ket = '>';
    else{
      std::cerr << "could not determine brackets" << std::endl;
      throw std::runtime_error("Bracket not recognised");
    }
    
    int level = 1;
    for( ; it < s.length() && level != 0; ++it){
      if (s[it] == bra) ++level;
      else if (s[it] == ket) --level;
    }
    
    if (level) return std::string::npos;
    else return it;
  }

  /**
   * parse a range into an
   * index array
   */
  void parse_range(std::string s, std::vector<int> & range, int x)
  {
    std::map<std::string, int> var;
    var["x"] = x;
    
    ExpressionParser<int> ep;
    
    std::string::size_type it = s.find(',');

    std::string rest;
    if (it != std::string::npos){
      rest = s.substr(it+1, std::string::npos);
      s = s.substr(0, it);
    }
    else{
      rest = "";
    }

    // parse the single part
    it = s.find('-');
    if (it != std::string::npos){
      int sr, er;
      sr = ep.parse_expression(s.substr(0, it), var);
      er = ep.parse_expression(s.substr(it+1, std::string::npos), var);
      for(int i=sr; i<=er; ++i)
	range.push_back(i);
    }
    else{
      int i = ep.parse_expression(s, var);
      range.push_back(i);
    }

    if (rest != "")
      parse_range(rest, range, x);
  }
}
