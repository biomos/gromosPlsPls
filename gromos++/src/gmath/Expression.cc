// gmath_Expression.cc

#include "Expression.h"
#include "../gromos/Exception.h"
#include <vector>
#include <sstream>
#include <iostream>

namespace gmath
{
  Expression::Expression(string s)
  {
    this->setExpression(s);
  } 

  Expression::Expression(string s, vector<double>& v)
  {
    this->setExpression(s);
    this->setValues(v);
  }

  void Expression::setExpression(string s)
  {
    vector<string> tokens;
    Expression::Tokenize(s, tokens, " ");
    d_val.resize(tokens.size());
    d_op.resize(tokens.size());
  
    for(unsigned int i=0; i< tokens.size(); i++){
      if(!allowed_token(tokens[i])){
	d_val[i]=atof(tokens[i].c_str());
	if(d_val[i]==0 && tokens[i]!="0"){
	  ostringstream os;
	  os << "Parse error: Do not know how to treat '" << tokens[i] 
	     << "' in expression. Allowed characters are +, -, *, /, (, )"
	     << ", a variable a<number>, the operators cos, sin and exp "
	     << "and any numbers" << endl;
	  throw(gromos::Exception("Expression", os.str()));
	}
	for(unsigned int l=1; l< tokens[i].size(); l++){
	  if(allowed_token(tokens[i].substr(l,l+1))){
	    ostringstream os;
	    os << "Parse error: space expected in " << tokens[i] << endl;
	    throw(gromos::Exception("Expression", os.str()));
	  }
	}
      }
      d_op[i]=tokens[i];
    }
    d_new=1;
    
    // do some tests
    if(d_op.size()==0)
      throw(gromos::Exception("Error in Expression", 
			      "Expression empty.\n"));
    if(d_op.size()==1){
      if(!is_number(d_op[0])){
	ostringstream os;
	os << "An expression that consists of only one "
	   << "token, cannot be an operator or a function (" << d_op[0] 
	   << ")." << endl;
	throw(gromos::Exception("Error in Expression", os.str()));
      }
    }
    else{
      int count_open_brackets=0;
      int count_close_brackets=0;
      if(d_op[0]=="(") count_open_brackets++;
      if(d_op[1]==")") count_close_brackets++;
      for(unsigned int i=1; i< d_op.size(); i++){
	if(is_number(d_op[i])){
	  //two numbers after each other
	  if(is_number(d_op[i-1])){
	    ostringstream os;
	    os << "Two numbers (" << d_op[i-1] << " and "
	       << d_op[i] << ") should be divided by an operator." << endl;
	    throw(gromos::Exception("Error in Expression", os.str()));
	  }
	  //a number after a bracket
	  if(d_op[i-1]==")"){
	    ostringstream os;
	    os << "A closing bracket cannot be followed "
	       << "by a number (" << d_op[i] << "). Put an operator in "
	       << "between." << endl;
	    throw(gromos::Exception("Error in Expression", os.str()));
	  }
	}
	// two operators after each other
	if(is_operator(d_op[i]) && is_operator(d_op[i-1])){
	  ostringstream os;
	  os << "An operator (" << d_op[i-1] 
	     << ") cannot be followed by another operator (" << d_op[i]
	     << ")." << endl;
	  throw(gromos::Exception("Error in Expression", os.str()));
	}
        // a number before a function
	if(is_number(d_op[i-1]) && is_function(d_op[i])){
	  ostringstream os;
	  os << "A number (" << d_op[i-1]
	     << ") cannot be followed by a function (" << d_op[i]
	     << ")." << endl;
	  throw(gromos::Exception("Error in Expression", os.str()));
	}
        // a function should be followed by a number or an opening bracket
	if(is_function(d_op[i-1]))
	  if ( (!is_number(d_op[i])) && d_op[i]!="("){
	    ostringstream os;
	    os << "A function (" << d_op[i-1]
	       << ") should be followed by a number or an opening bracket, "
	       << "not by (" << d_op[i] << ")." << endl;
	    throw(gromos::Exception("Error in Expression", os.str()));
	  }
	// some things about brackets
	if(d_op[i]=="(") {
	  count_open_brackets++;
	  // Two brackets )( after each other
	  if(d_op[i-1]==")"){
	    ostringstream os;
	    os << "Closing and opening brackets cannot "
	       << "come after each other. Put an operator in between." <<endl;
	    throw(gromos::Exception("Error in Expression", os.str()));
	  }
	  // number followed by a bracket
	  if(is_number(d_op[i-1])){
	    ostringstream os;
	    os << "An opening bracket cannot be "
	       << "preceded by a number (" << d_op[i-1] << "). Put an "
	       << "operator in between." << endl;
	    throw(gromos::Exception("Error in Expression", os.str()));
	  }
	}
	if(d_op[i]==")") {
	  count_close_brackets++;
	  // Two brackets () directly after eacht other
	  if(d_op[i-1]=="("){
	    ostringstream os;
	    os << "Opening and closing brackets cannot "
	       << "come after each other. Put something in between." << endl;
	    throw(gromos::Exception("Error in Expression", os.str()));
	  }
	}
      }
      if(count_open_brackets!=count_close_brackets){
	ostringstream os;
	os << "Found " << count_open_brackets 
	   << " opening brackets and " << count_close_brackets 
	   << " closing brackets." << endl;
	throw(gromos::Exception("Error in Expression", os.str()));
      }
    }
  }

  void Expression::setValues(vector<double>& v)
  {
    // first check how many we expect
    int max=-1;
    vector<int> vars;
    vector<int> varNumber;
  
    for(unsigned int i=0; i< d_op.size(); i++){
      if(d_op[i].substr(0,1)=="a"){
        vars.push_back(i);
	int num=atoi(d_op[i].substr(1, d_op[i].size()-1).c_str())-1;
	varNumber.push_back(num);
	if(num>max) max=num;
      }
    }
    if(max+1!=int(v.size())){
      ostringstream os;
      os << "Incorrect number of values supplied in function setValues.\n"
	 << "Expected " << max+1 << " values. Got " << v.size() << endl;
      throw(gromos::Exception("Expression", os.str()));
    }

    for(unsigned int i=0; i< vars.size(); i++){
      d_val[vars[i]]=v[varNumber[i]];
    }
    d_new=1;
  }

  
  void Expression::writeExpression(ostream& os)
  {
    for(unsigned int i=0; i< d_op.size(); i++){
      os << d_op[i] << " ";
    }
    os << endl;
  }

  void Expression::writeExpressionValue(ostream& os)
  {
    for(unsigned int i=0; i< d_op.size(); i++){
      if(is_number(d_op[i]))
	os << d_val[i] << " ";
      else
	os << d_op[i] << " ";
    }
    os << endl;
  }

  double Expression::value()
  {
    if(!d_new)
      return d_result;
    else{
      // as the calculation process changes the vectors with values
      // and operators, we have to make copies of these 
      vector<string> op;
      vector<double> val;
      for(unsigned int i=0; i< d_op.size(); i++){
	op.push_back(d_op[i]);
	val.push_back(d_val[i]);
      }
    
      d_result=calc(op, val, 0,d_op.size());
    
      d_new=0;
      
      return d_result;
    }
  }
  
  void Expression::Tokenize(const string& str,
			    vector<string>& tokens,
			    const string& delimiters)
  {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
  }

  double Expression::calc(std::vector<string>& op, std::vector<double>& val, 
			  int f, int t)
  {
    int i=f;
    int option=-1;

    // do we have brackets?
    // first get rid of all brackets
    int brackets=1;
    while(brackets){
      int br_open_first=f;
      int br_close_last=t;
    
      Expression::find_bracket(op, br_open_first, br_close_last);
      int tobesubstracted=br_close_last-br_open_first;
      if(br_open_first!=-1){
	double c=Expression::calc(op, val, br_open_first+1,br_close_last);
	br_open_first=f;
	br_close_last=val.size()-1;
	Expression::find_bracket(op, br_open_first, br_close_last);

	// now remove the elements between the brackets
	vector<double>::iterator begin= val.begin() + br_open_first,
	  end = val.begin()+br_close_last;
	vector<string>::iterator op_beg= op.begin() + br_open_first,
	  op_end = op.begin()+br_close_last;
      
	val.erase(begin, end);
	op.erase(op_beg,  op_end);
	val[br_open_first]=c;
	op[br_open_first]=" ";
	t-=tobesubstracted;
      }
      else brackets=0;
    }
    
    // special function
    special_functions(op, val, f, t);

    // now do the calculation, search for an operator
    // we prefer a + or a -, but keep the position of * or / in case we do
    // not find the first    
    while( i<t && op[i]!="+" && op[i]!="-" ) {
      if(op[i]=="*" || op[i]=="/") option = i;
      i++;
    }
    if(i==t) {
      if(option==-1){
	// No operator, so just return the value
        return val[i-1];
      }
      else i=option;
    }

    double a=calc(op, val, f, i);
    double b=calc(op, val, i+1, t);
    if(op[i]=="*")
	return a*b;
    if(op[i]=="/")
	return a/b;
    if(op[i]=="-")
	return a-b;
    if(op[i]=="+")
	return a+b;
    return 0;
  }
  bool Expression::allowed_token(string s)
  {
      return s=="+" || s=="-" || s=="*" || s=="/" || s=="(" ||s==")" ||
	  s[0]=='a' || s=="cos" || s=="sin" || s=="exp";
  }

  bool Expression::is_operator(string s)
  {
    return s=="+" || s=="-" || s=="*" || s=="/";
  }
  bool Expression::is_variable(string s)
  {
    return s[0]=='a';
  }
  bool Expression::is_function(string s, double on, double& res)
  {
    if(s=="cos"){ res=cos(on); return true;}
    if(s=="sin"){ res=sin(on); return true;}
    if(s=="exp"){ res=exp(on); return true;}
    return false;
  }
  bool Expression::is_function(string s)
  {
    return s=="cos" || s=="sin" || s=="exp";
  }
  bool Expression::is_bracket(string s)
  {
    return s=="(" || s==")";
  }
  bool Expression::is_number(string s)
  {
      return !is_operator(s) && !is_bracket(s) && !is_function(s);
  }
  void Expression::find_bracket(vector<string>& op, int &first, int &last)
  {
    int f=first;
    int t=last;
    int br_open_cnt=0;
    for(int j=f; j<t; j++){
      if(op[j]=="(") {
        if(!br_open_cnt) first=j;
	br_open_cnt++;
      }
    }
    if(br_open_cnt){
      // find the matching bracket to the last one
      int countbracket=1, k=first+1;
      while(countbracket!=0 && k<=last){
        if(op[k]=="(") countbracket++;
        if(op[k]==")") countbracket--;
	k++;
      }
      last=k-1;
    }
    else{
      first=-1;
      last=-1;
    }
  }

  void Expression::special_functions(vector<string>& op, vector<double>& val, 
			 int f, int& t)
  {
    double c=0;
    for(int j=f; j<t; j++){
      if(is_function(op[j], val[j+1], c)){
	// now remove the element j+1 
	vector<double>::iterator begin= val.begin() + j,
	    end = val.begin()+j+1;
	vector<string>::iterator op_beg= op.begin() + j,
	    op_end = op.begin()+j+1;
	  
	  
	val.erase(begin, end);
	op.erase(op_beg,  op_end);
	val[j]=c;
	op[j]=" ";
	t--;
      }
    }
  }
}
  

