// ExpressionParser.h

#ifndef INCLUDED_UTILS_EXPRESSIONPARSER
#define INCLUDED_UTILS_EXPRESSIONPARSER

namespace utils
{

  //////////////////////////////////////////////////
  // expression parser
  //////////////////////////////////////////////////
  template<typename T>
  class ExpressionParser
  {
  public:
    enum operation_enum{
      // unary expressions
      op_sin = 0,
      op_exp = 1,
      op_ln = 2,
      // binary expressions
      op_mul = 3,
      op_div = 4,
      op_add = 5,
      op_sub = 6,
      op_undef = 100
    };

    /**
     * Constructor
     */
    ExpressionParser();
    
    /**
     * parse an expression
     */
    T parse_expression(std::string s,
		       std::map<std::string, T> & var);

    bool parse_unary(std::string & s, 
		     std::map<std::string, T> & var,
		     T & res);
    
  private:
    T _parse_expression(T arg1, operation_enum op1,
			std::string rest,
			std::map<std::string, T> & var);
    
    T do_operation(T arg1, T arg2, operation_enum op);

    std::map<std::string, operation_enum> d_op;
    const std::string d_op_string;
  };


  
}

#include "ExpressionParser.cc"

#endif
