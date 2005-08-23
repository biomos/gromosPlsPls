// ExpressionParser.cc

#include "parse.h"
#include <cmath>

namespace utils
{

  inline int sin(int x)
  {
    return int(::sin(double(x)));
  }
  inline int exp(int x)
  {
    return int(::exp(double(x)));
  }
  
  

  template<typename T>
  ExpressionParser<T>::ExpressionParser()
    : d_op_string("*/+-")
  {
    d_op["sin"] = op_sin;
    d_op["exp"] = op_exp;
    d_op["ln"] = op_ln;
    d_op["*"] = op_mul;
    d_op["/"] = op_div;
    d_op["+"] = op_add;
    d_op["-"] = op_sub;
  }
  
  template<typename T>
  T ExpressionParser<T>::parse_expression(std::string s,
					  std::map<std::string, T> & var)
  {

    T arg1;
    operation_enum op1 = op_undef;
    std::string rest;
    
    bool unary = parse_unary(s, var, arg1);
        
    std::string::size_type op1_it = s.find_first_of(d_op_string);
    if (op1_it != std::string::npos){
      
      std::string sop1 = s.substr(op1_it, 1);
      if (d_op.count(sop1) <= 0)
	throw std::runtime_error("unsupported operation!");
      
      op1 = d_op[sop1];
      rest = s.substr(op1_it+1, std::string::npos);
      s = s.substr(0, op1_it);
    }
    
    if (!unary){
      std::string sarg1;
      std::istringstream is(s);
      is >> sarg1;
      if (var.count(sarg1)){
	arg1 = var[sarg1];
      }
      else{
	is.clear();
	is.str(s);
	if (!(is >> arg1))
	  throw std::runtime_error("could not parse number");
      }
    }

    if (op1_it == std::string::npos) return arg1;

    // now we have arg1, op1 and rest
    return _parse_expression(arg1, op1, rest, var);
  }
  
  template<typename T>
  T ExpressionParser<T>::_parse_expression(T arg1, operation_enum op1,
					   std::string s,
					   std::map<std::string, T> & var)
  {
    // lookup next argument / operator pair
    T arg2;
    operation_enum op2 = op_undef;
    std::string rest;
    
    bool unary = parse_unary(s, var, arg2);
    
    std::string::size_type op2_it = s.find_first_of(d_op_string);
    if (op2_it != std::string::npos){
      
      std::string sop2 = s.substr(op2_it, 1);
      if (d_op.count(sop2) <= 0)
	throw std::runtime_error("unsupported operation!");
      
      op2 = d_op[sop2];
      rest = s.substr(op2_it+1, std::string::npos);
      s = s.substr(0, op2_it);
    }
    
    if (!unary){
      std::string sarg2;
      std::istringstream is(s);
      is >> sarg2;
      if (var.count(sarg2)){
	arg2 = var[sarg2];
      }
      else{
	is.clear();
	is.str(s);
	if (!(is >> arg2))
	  throw std::runtime_error("could not parse number");
      }
    }

    if (op2_it == std::string::npos) return do_operation(arg1, arg2, op1);
    
    // now we have arg1 op1 arg2 op2 and rest
    if (op1 < op2){
      return _parse_expression(do_operation(arg1, arg2, op1), op2, rest, var);
    }
    else{
      return do_operation(arg1, _parse_expression(arg2, op2, rest, var), op1);
    }
  }

  template<typename T>
  T ExpressionParser<T>::do_operation(T arg1, T arg2, operation_enum op)
  {
    switch(op){
      case op_mul : return arg1 * arg2;
      case op_div : return arg1 / arg2;
      case op_add : return arg1 + arg2;
      case op_sub : return arg1 - arg2;
      default:
	throw std::runtime_error("unsupported operation");
    }
  }

  template<typename T>
  bool ExpressionParser<T>::parse_unary(std::string & s,
					std::map<std::string, T> & var,
					T & res)
  {
    // check whether expression starts with bracket
    std::string::size_type bra = s.find_first_not_of(" ");

    if (s[bra] == '('){
      std::string::size_type ket = find_matching_bracket(s, '(', bra+1);
      if (ket == std::string::npos){
	throw std::runtime_error("could not find matching bracket");
      }
      
      res = parse_expression(s.substr(bra+1, ket-bra-2), var);

      if (s.length() > ket+1)
	s = s.substr(ket+1, std::string::npos);
      else s = "";
      
      return true;
    }
    else if (s.substr(bra, 4) == "sin("){
      std::string::size_type ket = find_matching_bracket(s, '(', bra+3);
      res = sin(parse_expression(s.substr(bra + 4, ket - bra - 5), var));
      
      s = s.substr(ket+1, std::string::npos);
      return true;
    }
    else if (s.substr(bra, 4) == "exp("){
      std::string::size_type ket = find_matching_bracket(s, '(', bra+3);
      res = exp(parse_expression(s.substr(bra + 4, ket - bra - 5), var));
      
      s = s.substr(ket+1, std::string::npos);
      return true;
    }
    else if (s.substr(bra, 4) == "ln("){
      std::string::size_type ket = find_matching_bracket(s, '(', bra+2);
      res = exp(parse_expression(s.substr(bra + 3, ket - bra - 4), var));
      
      s = s.substr(ket+1, std::string::npos);
      return true;
    }

    return false;
  }
  
}
