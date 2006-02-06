// ExpressionParser.h

#ifndef INCLUDED_UTILS_EXPRESSIONPARSER
#define INCLUDED_UTILS_EXPRESSIONPARSER

#include <stack>
#include <map>
#include <stdexcept>

#include "../bound/Boundary.h"
#include "Value.h"

namespace utils
{
  enum operation_enum{
    // functions
    // - basic
    op_abs = 10,
    // - trigonometric:
    op_sin = 20,
    op_cos = 21,
    op_tan = 22,
    op_asin = 23,
    op_acos = 24,
    op_atan = 25,
    // - transcendent:
    op_exp = 30,
    op_ln = 31,
    // - vector:
    op_abs2 = 51,
    op_dot = 52,
    op_cross = 53,
    op_ni = 54,
    // unary expressions
    op_unary = 100,
    op_umin = 101,
    op_uplus = 102,
    op_not = 103,
    // binary expressions
    op_binary = 200,
    op_mul = 201,
    op_div = 202,
    op_add = 203,
    op_sub = 204,
    // logical expressions
    op_logical = 300,
    op_eq = 301,
    op_neq = 302,
    op_less = 303,
    op_greater = 304,
    op_lesseq = 305,
    op_greatereq = 306,
    op_and = 400,
    op_or = 401,
    // ternary operator
    op_ternary = 500,
    op_condition = 501,
    // no op
    op_undef = 1000
  };

  template<typename T>
  class ExpressionParser;

  template<typename T>
  class ValueTraits;
  
  template<>
  class ValueTraits<int>;
  
  template<>
  class ValueTraits<Value>;

  //////////////////////////////////////////////////
  // expression parser
  //////////////////////////////////////////////////
  template<typename T>
  class ExpressionParser
  {
  public:

    enum expr_enum{
      expr_value,
      expr_operator,
      expr_function,
      expr_variable
    };
    
    struct expr_struct
    {
      expr_struct(T value) 
	: type(expr_value), value(value), op(op_undef), name("") {}
      expr_struct(operation_enum op) 
	: type(expr_operator), value(T()), op(op), name("") 
      {
	if (op < op_unary) type = expr_function;
      }
      expr_struct(std::string s) 
	: type(expr_variable), value(T()), op(op_undef), name(s) {}
      

      expr_enum type;
      T value;
      operation_enum op;
      std::string name;
      
      std::string toString()
      {
	std::ostringstream os;
	
	switch(type){
	  case expr_value: os <<    "value: " << value; break;
	  case expr_operator: os << "op:    " << op; break;
	  case expr_function: os << "f:     " << op; break;
	  case expr_variable: os << "var:   " << name; break;
	  default: os << "-- unknown --";
	}
	return os.str();
      }
    };

    /**
     * Constructor
     */
    ExpressionParser();

    ExpressionParser(gcore::System * sys,
		     bound::Boundary * pbc);
    
    void init_op();
    
    /**
     * parse an expression
     */
    T parse_expression(std::string s,
		       std::map<std::string, T> & var);

    void parse_expression(std::string s,
			  std::map<std::string, T> & var,
			  std::vector<expr_struct> & expr);

    void print_expression(std::vector<expr_struct> & expr,
			  std::ostream & os = std::cout);

    bool parse_unary(std::string & s, 
		     std::map<std::string, T> & var,
		     T & res);
    
    T calculate(std::vector<expr_struct> & expr,
		std::map<std::string, T> & var);

    void calculate(std::string name,
		   std::map<std::string, std::vector<expr_struct> > & expr,
		   std::map<std::string, T> & var);

    void calculate(std::map<std::string, std::vector<expr_struct> > & expr,
		   std::map<std::string, T> & var);

    void do_general_function(operation_enum op);
    void do_trigonometric_function(operation_enum op);
    void do_transcendent_function(operation_enum op);
    void do_vector_function(operation_enum op);
    void do_logical_operation(operation_enum op);

  private:
    T _parse_expression(T arg1, operation_enum op1,
			std::string rest,
			std::map<std::string, T> & var);

    void _parse_token(operation_enum op,
		      std::string s,
		      std::string::size_type & it,
		      std::map<std::string, T> & var,
		      std::vector<expr_struct> & expr);

    bool _parse_function(std::string s,
			 std::string::size_type & it,
			 std::map<std::string, T> & var,
			 std::vector<expr_struct> & expr);

    void _parse_value(std::string s,
		     std::string::size_type & it,
		     std::map<std::string, T> & var,
		     std::vector<expr_struct> & expr);
    
    void _commit_value(std::string s,
		      std::string::size_type & it,
		      std::map<std::string, T> & var,
		      std::vector<expr_struct> & expr);
    
    void _commit_operator(operation_enum op,
			  std::vector<expr_struct> & expr);
    
    operation_enum _parse_unary_operator(std::string s,
					 std::string::size_type & it);

    operation_enum _parse_operator(std::string s,
				   std::string::size_type & it);

    T do_operation(operation_enum op);
    void do_function(operation_enum op);

    T do_operation(T arg1, T arg2, operation_enum op);

    std::map<std::string, operation_enum> d_op;
    std::string d_op_string;

    std::stack<T> d_stack;

    ValueTraits<T> d_value_traits;
  };

  
  //////////////////////////////////////////////////
  // traits class for values
  //////////////////////////////////////////////////
  template<typename T>
  class ValueTraits
  {
  public:
    ValueTraits(gcore::System * sys,
		bound::Boundary *pbc) : d_sys(sys), d_pbc(pbc) {}
    
    ValueTraits() : d_sys(NULL), d_pbc(NULL) {}

    T parseValue(std::string s, std::map<std::string, T> & var)
    {
      std::istringstream is(s);
      T t;
      if (!(is >> t))
	throw std::runtime_error("could not read value (" + s + ")");
      return t;
    }
    
    static void do_function(operation_enum op,
			    ExpressionParser<T> & ep)
    {
      ep.do_general_function(op);
      ep.do_trigonometric_function(op);
      ep.do_transcendent_function(op);
      // ep.do_vector_function(op);
    }

    static void do_operation(operation_enum op,
			     ExpressionParser<T> & ep)
    {
      ep.do_logical_operation(op);
    }
    

  private:
    gcore::System * d_sys;
    bound::Boundary * d_pbc;
  };

  // integer specialisation
  template<>
  class ValueTraits<int>
  {
  public:
    ValueTraits() {}

    int parseValue(std::string s, std::map<std::string, int> & var)
    {
      std::istringstream is(s);
      int t;
      if (!(is >> t))
	throw std::runtime_error("could not read value (" + s + ")");
      return t;
    }
    
    static void do_function(operation_enum op,
			    ExpressionParser<int> & ep)
    {
      ep.do_general_function(op);
    }

    static void do_operation(operation_enum op,
			     ExpressionParser<int> & ep)
    {
      ep.do_logical_operation(op);
    }

  };

  // and finally the one for Values!
  template<>
  class ValueTraits<Value>
  {
  public:
    ValueTraits(gcore::System *sys,
		bound::Boundary *pbc) : d_sys(sys), d_pbc(pbc) {}
    
    // ValueTraits() : d_sys(NULL), d_pbc(NULL) {}

    Value parseValue(std::string s, std::map<std::string, Value> & var)
    {
      Value v;
      v.parse(s, var, *d_sys, d_pbc);
      return v;
    }
    
    static void do_function(operation_enum op,
			    ExpressionParser<Value> & ep)
    {
      ep.do_general_function(op);
      ep.do_trigonometric_function(op);
      ep.do_transcendent_function(op);
      ep.do_vector_function(op);
    }

    static void do_operation(operation_enum op,
			     ExpressionParser<Value> & ep)
    {
      ep.do_logical_operation(op);
    }

    bound::Boundary * pbc()
    {
      if (d_pbc == NULL) throw std::runtime_error("pbc is NULL");
      return d_pbc;
    }

    gcore::System * sys()
    {
      if (d_sys == NULL) throw std::runtime_error("sys is NULL");
      return d_sys;
    }

  private:
    gcore::System * d_sys;
    bound::Boundary * d_pbc;
  };

}

#include "ExpressionParser.cc"

#endif
