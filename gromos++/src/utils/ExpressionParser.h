// ExpressionParser.h

#ifndef INCLUDED_UTILS_EXPRESSIONPARSER
#define INCLUDED_UTILS_EXPRESSIONPARSER

#include <stack>

namespace utils
{
  enum operation_enum{
    // functions
    op_sin = 0,
    op_exp = 1,
    op_ln = 2,
    // unary expressions
    op_unary = 100,
    op_umin = 101,
    op_uplus = 102,
    // binary expressions
    op_binary = 200,
    op_mul = 201,
    op_div = 202,
    op_add = 203,
    op_sub = 204,
    // no op
    op_undef = 1000
  };

  template<typename T>
  class ExpressionParser;

  template<typename T>
  class ValueTraits;
  
  template<>
  class ValueTraits<int>;
  
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

    ExpressionParser(gcore::System & sys,
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

    bool parse_unary(std::string & s, 
		     std::map<std::string, T> & var,
		     T & res);
    
    T calculate(std::vector<expr_struct> & expr,
		std::map<std::string, T> & var);


    void do_trigonometric_function(operation_enum op);
    void do_transcendent_function(operation_enum op);
    void do_vector_function(operation_enum op);

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

    int _commit_value(std::string s,
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

    T parseValue(std::string s)
    {
      std::istringstream is(s);
      T t;
      if (!(is >> t))
	throw std::runtime_error("could not read value");
      return t;
    }
    
    static void do_function(operation_enum op,
			    ExpressionParser<T> & ep)
    {
      ep.do_trigonometric_function(op);
      ep.do_transcendent_function(op);
      // ep.do_vector_function(op);
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

    static const bool expr_vector_functions = false;
    static const bool expr_transcendent_functions = false;
    static const bool expr_trigonometric_functions = false;

    int parseValue(std::string s)
    {
      std::istringstream is(s);
      int t;
      if (!(is >> t))
	throw std::runtime_error("could not read value");
      return t;
    }
    
    static void do_function(operation_enum op,
			    ExpressionParser<int> & ep)
    {
    }

  };
}

#include "ExpressionParser.cc"

#endif
