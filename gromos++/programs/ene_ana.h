void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters = " ");

void print(gmath::stat &p, string s, double time, double dt);

class Expression
{
  vector<double> d_val;
  vector<string> d_op;
  int d_new;
  double d_result;

 public:
  Expression(string s);
  Expression(string s, vector<double>& v);
  
  void setExpression(string s);
  void setValues(vector<double>& v);
  double value();
  void writeExpression(ostream& os = std::cout);
  void writeExpressionValue(ostream& os = std::cout);
  
private:
  void Tokenize(const string& str,
		vector<string>& tokens,
		const string& delimiters = " ");
  double calc(vector<string>& op, vector<double>& val, int f, int t);
  void find_bracket(vector<string>&op, int &first, int &last);
};

class energy_trajectory
{
  vector<double> d_data;
  vector<Expression> d_e;
  vector<vector <int> > d_dep;
  vector<double> d_calc;
  vector<bool> d_recalc;
  
  map<string, int> d_map;
  int d_frame;
  int d_negr;
public:
  energy_trajectory();
  double operator[](int i);
  double operator[](string s);
  void read_frame(Ginstream& is);
  void addKnown(string s, string v);
  int index(string s);
  
  void write_map(ostream& os=cout);
 private:
  string back_index(int i);
  
  
  
};
typedef map<string, int>::value_type MP;
const int num_ener=22;
const int num_eneres=28;
const int num_volprt=48;
const int unknownvariable=-1000;


void set_standards(energy_trajectory &e, double mass);
void read_library(string name, energy_trajectory& e);

/**********************************************************
 implementation
**********************************************************/
void print(gmath::stat &p, string s, double time, double dt)
{
  ostrstream os;
  os << s << ".out" << ends;
  ofstream fout(os.str());
  fout << "#"
       << setw(9) << "time"
       << setw(14) << s
       << endl;
  for(int i=0; i< p.n(); i++){
    fout << setw(10) << time+i*dt
	 << setw(14) << p.val(i)
	 << endl;
  }
  fout.close();
  // and print the averages etc to cout
  cout << setw(10) << s
       << setw(14) << p.ave()
       << setw(14) << p.rmsd()
       << setw(14) << p.ee()
       << endl;
}

energy_trajectory::energy_trajectory()
{
  d_frame=0;
  int i=0;
  // set the really standard names: ENER, ENERES, VOLPRT  
  for(; i<num_ener; i++){
    ostrstream os;
    os << "ENER[" << i+1 << "]" << ends;
    d_map.insert(MP(os.str(), i));
  }
  for(; i<num_eneres; i++){
    ostrstream os;
    os << "ENERES[" << i-num_ener+1 << "]" << ends;
    d_map.insert(MP(os.str(), i));
  }
  for(; i<num_volprt; i++){
    ostrstream os;
    os << "VOLPRT[" << i-num_eneres+1 << "]" << ends;
    d_map.insert(MP(os.str(), i));
  }
}

int energy_trajectory::index(string s)
{
  map<string, int>::const_iterator iter=d_map.find(s);
  if(iter!=d_map.end()) return iter->second;
  
  // otherwise it might be one of the charge groups
  string sub=s.substr(0,6);
  int g_num = atoi(s.substr(s.find("[")+1, s.find("]")).c_str())-1;

  if(sub=="ENERLJ") return num_volprt + 4*g_num;
  if(sub=="ENERCL") return num_volprt + 4*g_num + 1;
  if(sub=="ENERRF") return num_volprt + 4*g_num + 2;
  if(sub=="ENERRC") return num_volprt + 4*g_num + 3;

  return unknownvariable;
}

    
double energy_trajectory::operator[](int i)
{
  if(i!=unknownvariable){
    
    if(i>=0 && i<num_volprt+d_negr)
      return d_data[i];
    else if (i<0){
      
      //calculate
      int ind= -i-1;
      if(d_recalc[ind]){
	
	vector<double> v;
	for(unsigned int k=0; k< d_dep[ind].size(); k++)
	  v.push_back((*this)[d_dep[ind][k]]);
	d_e[ind].setValues(v);
	// d_e[ind].writeExpression(cout);
	// d_e[ind].writeExpressionValue(cout);
	d_recalc[ind]=false;
	d_calc[ind]=d_e[ind].value();
      }
      return d_calc[ind];
    }
  }
  
  cerr << "Trying to access an unkwown element: " << i << endl;
  
  return 0.0;
}
double energy_trajectory::operator[](string s)
{
  return (*this)[this->index(s)];
}


void energy_trajectory::addKnown(string s, string e)
{
  // do we know this s already
  int i=this->index(s);
  int known=0;
  if(i != unknownvariable) known = 1;

  vector<string> v;
  Tokenize(e, v);

  int ind=this->index(v[0]);
  
  // if it is just a one-to-one mapping, add it to the map (or reset the map)
  if(v.size()==1 && ind != unknownvariable){
    if(!known)
      d_map.insert(MP(s, ind) );
    else
      d_map[s]=ind;
  }
  else{
    // we have an expression
    ostrstream os;
    int varcount=1;
    vector<int> dep;
    
    for(unsigned int j=0; j< v.size(); j++){
      int var=this->index(v[j]);
      // check if it is a known variable
      if(var!=unknownvariable) {
	os << "a" << varcount;
	dep.push_back(var);
	varcount++;
      }
      // otherwise it must be either a number or an operator
      else os << v[j];
      os << " ";
    }
    os << ends;
    
    Expression e(os.str());
    
    if(!known){
      d_e.push_back(e);
      d_dep.push_back(dep);
      d_calc.push_back(0.0);
      d_recalc.push_back(true);
      d_map.insert(MP(s, -1*d_e.size()));
    }
    else{
      d_e[i].setExpression(os.str());
      d_dep[i].resize(0);
      for(unsigned int k=0; k<dep.size(); k++)
	d_dep[i].push_back(dep[k]);
    }
  }
}

    
void energy_trajectory::read_frame(Ginstream& is)
{
  if(!d_frame)
    d_data.resize(num_volprt,0.0);
  // some temporary things to store the energygroup energies
  vector<double> energrp;
  string sdum;
  while(sdum!="ENERGY") { is >> sdum; }
  
  
  //first read in the energy block
  int i=0;
  for(; i<num_ener; i++)
    is >> d_data[i];
  //the eneres block
  for(; i<num_eneres; i++)
    is >> d_data[i];
  // read the number of energy groups
  int negr, size;
  is >> negr;
  size=4*negr*(negr+1)/2;
  if(!d_frame){
    d_negr=size;
  }
  else if(size!=d_negr)
    cerr << "\nNumber of energy groups is not constant!\n";
  energrp.resize(size, 0.0);  	
  
  for(int j=0; j<size; j++){
    is >> energrp[j];
  }
  while(sdum!="VOLUMEPRESSURE") {is >> sdum; }
  
  
  // now read in the volume pressure block
  for(;i<num_volprt; i++)
    is >> d_data[i];
  // possibly resize the d_data array
  if(!d_frame)
    d_data.resize(num_volprt+4*size);
  for(int j=0; j<size; j++){
    d_data[num_volprt+j]=energrp[j];
  }
  d_frame++;
  is >> sdum >> sdum;

  //set all things that need to be calculated to be uncalculated
  for(unsigned int i=0; i<d_recalc.size(); i++) d_recalc[i]=true;
  
}

string energy_trajectory::back_index(int i)
{
  map<string, int>::const_iterator iter=d_map.begin(), to=d_map.end();
  for(; iter!=to; ++iter)
    if(iter->second == i) return iter->first;
  if(i>=num_volprt){
    ostrstream os;
    
    int j=i-num_volprt;
    switch(j%4){
      case 0: os << "ENERLJ[" << int(j/4)+1 << "]" << ends;
	break;
      case 1: os << "ENERCL[" << int(j/4)+1 << "]" << ends;
	break;
      case 2: os << "ENERRF[" << int(j/4)+1 << "]" << ends;
	break;
      case 3: os << "ENERRC[" << int(j/4)+1 << "]" << ends;
	break;
    }
    return os.str();
  }
  
  return "unknown";
}

void energy_trajectory::write_map(ostream& os)
{
  map<string, int>::const_iterator iter=d_map.begin(), to=d_map.end();
  for(;iter!=to; ++iter){
    os << setw(12) << iter->first << " = ";
    os << "data[" << iter->second << "]";
    string nm=this->back_index(iter->second);
    if(iter->first != nm)
      os << "  (" << nm << ")"; 
    os << endl;
    if(iter->second < 0){
      os << setw(15) << "= ";
      int i=-1-iter->second;
      d_e[i].writeExpression(os);
      if(d_dep[i].size()!=0){
	os << setw(20) << "with:" << endl;
	for(unsigned int j=0; j<d_dep[i].size(); j++)
	  os << setw(20) << "a" << j+1 << " = data[" 
	     << d_dep[i][j] << "]  (" << this->back_index(d_dep[i][j])
	     << ")" << endl;
      }
    }
  }
}

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
  Expression::Tokenize(s, tokens);
  string ops="+-*/()";
  d_val.resize(tokens.size());
  d_op.resize(tokens.size());
  
  for(unsigned int i=0; i< tokens.size(); i++){
    if(ops.find(tokens[i], 0) == string::npos && tokens[i][0]!='a'){
      d_val[i]=atof(tokens[i].c_str());
      if(d_val[i]==0 && tokens[i]!="0")
	cerr << "Parse error: Do not know how to treat '" << tokens[i] 
	     << "' in expression. Allowed characters are " << ops
	     << ", a variable a<number> and any numbers" << endl;
      for(unsigned int l=0; l< ops.size(); l++){
	if(tokens[i].find(ops.substr(l,l+1),0) !=string::npos)
	  cerr << "Parse error: " << tokens[i] << " will be interpreted as "
	       << d_val[i] << endl;
      }
    }
    d_op[i]=tokens[i];
  }
  d_new=1;
  
  // do some tests
  if(d_op.size()==0)
    cerr<< "Error in expression: expression empty." << endl;
  if(d_op.size()==1){
    if(ops.find(d_op[0])!=string::npos)
      cerr << "Error in expression: An expression that consists of only one "
	   << "token, cannot be an operator (" << d_op[0] << ")." << endl;
  }
  else{
    string ops_red="+-*/()";
    string ops_ops="+-*/";
    int count_open_brackets=0;
    int count_close_brackets=0;
    if(d_op[0]=="(") count_open_brackets++;
    if(d_op[1]==")") count_close_brackets++;
    
    for(unsigned int i=1; i< d_op.size(); i++){
      if((ops_red.find(d_op[i].substr(0,1), 0) == string::npos)){
	//two numbers after each other
	if(ops_red.find(d_op[i-1].substr(0,1), 0) == string::npos)
	  cerr << "Error in expression: Two numbers (" << d_op[i-1] << " and "
	       << d_op[i] << ") should be divided by an operator." << endl;
	//a number after a bracket
	if(d_op[i-1]==")")
	  cerr << "Error in expression: A closing bracket cannot be followed "
	       << "by a number (" << d_op[i] << "). Put an operator in between."
	       << endl;
      }
      // two operators after each other
      if(ops_ops.find(d_op[i], 0) != string::npos && 
	 ops_ops.find(d_op[i-1], 0) != string::npos)
	cerr << "Error in expression: An operator (" << d_op[i-1] 
	     << ") cannot be followed by another operator (" << d_op[i]
	     << ")." << endl;
      // some things about brackets
      if(d_op[i]=="(") {
	count_open_brackets++;
	// Two brackets )( after each other
	if(d_op[i-1]==")")
	  cerr << "Error in expression: Closing and opening brackets cannot "
	       << "come after each other. Put an operator in between." <<endl;
	// number followed by a bracket
	if(ops_red.find(d_op[i-1].substr(0,1), 0) == string::npos)
	  cerr << "Error in expression: An opening bracket cannot be "
	       << "preceded by a number (" << d_op[i-1] << "). Put an "
	       << "operator in between." << endl;
      }
      if(d_op[i]==")") {
	count_close_brackets++;
	// Two brackets () directly after eacht other
	if(d_op[i-1]=="(")
	  cerr << "Error in expression: Opening and closing brackets cannot "
	       << "come after each other. Put something in between." << endl;
      }
    }
    if (count_open_brackets!=count_close_brackets)
      cerr << "Error in expression: found " << count_open_brackets 
	   << " opening brackets and " << count_close_brackets 
	   << " closing brackets." << endl;
  }

}

void Expression::setValues(vector<double>& v)
{
  // first check how many we expect
  int max=-1;
  vector<int> vars;
  vector<int> varNumber;
  
  for(unsigned int i=0; i< d_op.size(); i++){
    if (d_op[i].substr(0,1)=="a"){
      vars.push_back(i);
      int num=atoi(d_op[i].substr(1, d_op[i].size()-1).c_str())-1;
      varNumber.push_back(num);
      if(num>max) max=num;
    }
  }
  if(max+1!=int(v.size()))
    cerr << "Incorrect number of values supplied in function setValues.\n"
	 << "Expected " << max+1 << " values. Got " << v.size() << endl;
    
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
  string ops="+-*/()";
    
   
  for(unsigned int i=0; i< d_op.size(); i++){
    if(ops.find(d_op[i].substr(0,1),0) == string::npos)
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
void Tokenize(const string& str,
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

double Expression::calc(vector<string>& op, vector<double>& val, int f, int t)
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

  // now do the calculation    
  while( i<t && op[i]!="+" && op[i]!="-" ) {
    if(op[i]=="*" || op[i]=="/") option = i;
    i++;
  }
  
  if(i==t) {
    if(option==-1){
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


void set_standards(energy_trajectory &e, double mass)
{
  {
    ostrstream os;
    os << mass << ends;
    e.addKnown("MASS", os.str());
  }
  {
    ostrstream os;
    os << BOLTZ << ends;
    e.addKnown("BOLTZ", os.str());
  }
  
  e.addKnown("totene", "ENER[1]");
  e.addKnown("totkin", "ENER[2]");
  e.addKnown("totpot", "ENER[9]");
  e.addKnown("pressu", "VOLPRT[12] * 16.388453");
  e.addKnown("boxvol", "VOLPRT[8]");
  e.addKnown("densit", "MASS * 1.66056 / VOLPRT[8]");
  
}

void read_library(string name, energy_trajectory& e)
{
  Ginstream gin(name);
  string sdum;
  vector<string> data;
  int count=0;
  
  // first read in everything to a vector of strings
  while(sdum!="VARIABLES") gin >> sdum;
  while(sdum !="END") {
    gin >> sdum;
    data.push_back(sdum);
    count++;
  }

  // now search for the first appearance of "="
  for(unsigned int i=0; i< data.size(); i++){
    if(data[i]=="="){

      // search for either the next appearance or the end
      unsigned int to=i+1;
      for(; to < data.size(); to++) if(data[to]=="=") break;

      // parse the expression part into an ostrstream
      ostrstream os;
      for(unsigned int j=i+1; j< to-1; j++) os << " " << data[j]; 
      os << ends;
      
      e.addKnown(data[i-1], os.str());
    }
  }
}

  

