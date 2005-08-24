// ExpressionParser.cc

#include "parse.h"
#include <cmath>
#include <iomanip>

namespace utils
{

  //////////////////////////////////////////////////
  // ValueTraits
  //////////////////////////////////////////////////

  /*
  template<typename T>
  void ValueTraits<T>::do_function
  (operation_enum op,
   ExpressionParser<T> & ep)
  {
    ep.do_trigonometric_function(op);
    ep.do_transcendent_function(op);
    // ep.do_vector_function(op);
  }
  */

  /*
  template<>
  void ValueTraits<int>::do_function
  (operation_enum op,
   ExpressionParser<int> & ep)
  {
    // ep.do_trigonometric_function(op);
    // ep.do_transcendent_function(op);
    // ep.do_vector_function(op);
  }
  */
  

  //////////////////////////////////////////////////

  template<typename T>
  ExpressionParser<T>::ExpressionParser()
  {
    init_op();
  }
  
  template<typename T>
  ExpressionParser<T>::ExpressionParser(gcore::System &sys,
					bound::Boundary * pbc)
    : d_value_traits(sys, pbc)
  {
    init_op();
  }

  template<typename T>
  void ExpressionParser<T>::init_op()
  {
    d_op_string = "*/+-";

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
    std::string::size_type it = 0;
    std::vector<expr_struct> expr;

    _parse_token(op_undef, s, it, var, expr);

    // for(unsigned int i=0; i<expr.size(); ++i){
    // std::cout << std::setw(5) << i << ":    " << expr[i].toString() << std::endl;
    // }

    T res = calculate(expr, var);

    // std::cout << "result (new) = " << res << "\n" << std::endl;
    return res;
  }

  template<typename T>
  void ExpressionParser<T>::_parse_token(operation_enum op,
					 std::string s,
					 std::string::size_type & it,
					 std::map<std::string, T> & var,
					 std::vector<expr_struct> & expr)
  {
    // function, unary operator, value
    _parse_value(s, it, var, expr);
    
    operation_enum op2 = op_undef;
    if (it != std::string::npos)
      op2 = _parse_operator(s, it);
    // else std::cerr << "no second operator: end" << std::endl;
    
    if (op <= op2){
      if (op != op_undef)
	_commit_operator(op, expr);
      if (op2 != op_undef)
	_parse_token(op2, s, it, var, expr);
      // else std::cerr << "end reached" << std::endl;
    }
    else
      {
	_parse_token(op2, s, it, var, expr);
	if (op != op_undef)_commit_operator(op, expr);
      }
  }

  template<typename T>
  bool ExpressionParser<T>::_parse_function
  (
   std::string s, 
   std::string::size_type & it,
   std::map<std::string, T> & var,
   std::vector<expr_struct> & expr
   )
  {
    // std::cerr << "parsing function " 
    // << s.substr(it, std::string::npos) << std::endl;
    
    std::string::size_type bra = s.find_first_not_of(" ", it);
    operation_enum fop = op_undef;
    
    if (s[bra] == '('){
      bra += 1;
    }
    else if (s.substr(bra, 4) == "sin("){
      bra += 4;
      fop = op_sin;
    }
    else if (s.substr(bra, 4) == "exp("){
      bra += 4;
      fop = op_exp;
    }
    else if (s.substr(bra, 4) == "ln("){
      bra += 3;
      fop = op_ln;
    }
    else{
      return false;
    }
    
    std::string::size_type ket = find_matching_bracket(s, '(', bra+1);
    if (ket == std::string::npos){
      throw std::runtime_error("could not find matching bracket");
    }

    // std::cerr << "function: " << fop << std::endl;

    std::string::size_type fit = 0;
    _parse_token(op_undef, s.substr(bra, ket-bra-1), fit, var, expr);

    if (fop != op_undef)
      _commit_operator(fop, expr);

    if (s.length() > ket)
      it = ket;
    else it = std::string::npos;
    
    return true;
  }
  
  template<typename T>
  operation_enum
  ExpressionParser<T>::_parse_unary_operator
  (
   std::string s, 
   std::string::size_type & it
   )
  {
    // std::cerr << "parsing unary: " << s.substr(it, std::string::npos)
    // << std::endl;

    std::string::size_type bra = s.find_first_not_of(" ", it);
    if (bra == std::string::npos) return op_undef;
  
    // std::cerr << "\tunary: " << s[bra] << std::endl;
  
    if (s[bra] == '-'){ ++it; return op_umin; }
    if (s[bra] == '+'){ ++it; return op_uplus; }   

    return op_undef;
  }

  template<typename T>
  operation_enum
  ExpressionParser<T>::_parse_operator
  (
   std::string s, 
   std::string::size_type & it
   )
  {
    // std::cerr << "parse_operator: " << s.substr(it, std::string::npos)
    // << std::endl;
    
    std::string::size_type bra = s.find_first_of(d_op_string, it);
    if (bra == std::string::npos) return op_undef;
    
    std::string ops = s.substr(bra,1);
    // std::cerr << "ops: " << ops << std::endl;
    
    if (d_op.count(ops) == 0)
      throw std::runtime_error("operator undefined!");

    ++it;
    return d_op[ops];
  }

  template<typename T>
  void ExpressionParser<T>::_parse_value
  (
   std::string s, 
   std::string::size_type & it,
   std::map<std::string, T> & var,
   std::vector<expr_struct> & expr
   )
  {
    if (_parse_function(s, it, var, expr))
      return;

    operation_enum u_op = op_undef;
    u_op = _parse_unary_operator(s, it);

    if (u_op != op_undef){
      _parse_value(s, it, var, expr);
      _commit_operator(u_op, expr);
      return;
    }

    _commit_value(s, it, var, expr);
  }
  
  template<typename T>
  void ExpressionParser<T>::_commit_value
  (
   std::string s, 
   std::string::size_type & it,
   std::map<std::string, T> & var,
   std::vector<expr_struct> & expr
   )
  {
    // TODO: check for bracket / function here!

    std::string::size_type vit = s.find_first_of(d_op_string, it);
    
    // std::cerr << "committing " << s.substr(it, vit - it) << std::endl;
    
    // check whether it's a variable
    std::istringstream is(s.substr(it, vit-it));
    std::string name;
    is >> name;
    if (var.count(name) > 0){
      expr_struct e(name);
      expr.push_back(e);
    }
    else{
      expr_struct e(d_value_traits.parseValue(s.substr(it, vit - it)));
      expr.push_back(e);
    }
    
    it = vit;
    if (it >=  s.length()) it = std::string::npos;
  }

  template<typename T>
  void ExpressionParser<T>::_commit_operator
  (
   operation_enum op,
   std::vector<expr_struct> & expr
   )
  {
    // std::cerr << "committing operator " << op << std::endl;
    expr.push_back(expr_struct(op));
  }

  template<typename T>
  T ExpressionParser<T>::calculate
  (
   std::vector<expr_struct> & expr,
   std::map<std::string, T> & var
   )
  {
    for(unsigned int i=0; i<expr.size(); ++i){
      switch(expr[i].type){
	case expr_value: d_stack.push(expr[i].value); break;
	case expr_variable: d_stack.push(var[expr[i].name]); break;
	case expr_function: do_function(expr[i].op); break;
	case expr_operator: do_operation(expr[i].op); break;
	default: throw std::runtime_error("wrong expression type");
      }
    }
    
    T res = d_stack.top();
    d_stack.pop();
    return res;
  }

  template<typename T>
  T ExpressionParser<T>::do_operation(operation_enum op)
  {
    T res;
    if (op < op_unary) throw std::runtime_error("operator is function");
    if (op < op_binary){
      T arg = d_stack.top();
      d_stack.pop();
      switch(op){
	case op_uplus: break;
	case op_umin: res = -arg; break;
	default: throw std::runtime_error("unknown unary operator");
      }
    }
    else if (op < op_undef){
      T arg2 = d_stack.top();
      d_stack.pop();
      T arg1 = d_stack.top();
      d_stack.pop();
      switch(op){
	case op_add: res = arg1 + arg2; break;
	case op_sub: res = arg1 - arg2; break;
	case op_mul: res = arg1 * arg2; break;
	case op_div: res = arg1 / arg2; break;
	default: throw std::runtime_error("unknown binary operator");
      }
    }
    else{
      throw std::runtime_error("unknown / undef operator");
    }
    
    d_stack.push(res);
    return res;
  }

  template<typename T>
  void ExpressionParser<T>::do_function(operation_enum op)
  {
    ValueTraits<T>::do_function(op, *this);
  }
  
  template<typename T>
  void ExpressionParser<T>::do_trigonometric_function(operation_enum op)
  {
    if (op == op_undef) return;
    
    T res;
    T arg = d_stack.top();
    
    switch(op){
      case op_sin: res = sin(arg); break;
      default: return;
    }
    
    op = op_undef;
    d_stack.pop();
    d_stack.push(res);
  }

  template<typename T>
  void ExpressionParser<T>::do_transcendent_function(operation_enum op)
  {
    if (op == op_undef) return;
    
    T res;
    T arg = d_stack.top();
    
    switch(op){
      case op_exp: res = exp(arg); break;
      case op_ln:  res = log(arg); break;
      default: return;
    }
    
    op = op_undef;
    d_stack.pop();
    d_stack.push(res);
  }

  template<typename T>
  void ExpressionParser<T>::do_vector_function(operation_enum op)
  {
    if (op == op_undef) return;
    
    T res;
    
    switch(op){
      default: return;
    }

    op = op_undef;
    d_stack.push(res);
  }
}
