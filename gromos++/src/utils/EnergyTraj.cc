#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

#include "EnergyTraj.h"
#include "../gmath/Expression.h"
#include "../gromos/Exception.h"

using utils::EnergyTraj;

EnergyTraj::EnergyTraj()
{
  std::vector<int> numbers;
  numbers.push_back(22);
  numbers.push_back(6);
  numbers.push_back(20);
  numbers.push_back(22);
  numbers.push_back(9);
  init(numbers);
}

EnergyTraj::EnergyTraj(std::vector<int> numbers)
{
  init(numbers);
}

void EnergyTraj::init(std::vector<int> numbers)
{
  num_ener   = numbers[0];
  num_eneres = num_ener   + numbers[1];
  num_volprt = num_eneres + numbers[2];
  num_rlam   = num_volprt + 1;
  num_fren   = num_rlam   + numbers[3];
  num_enerf  = numbers[4];
  unknownvariable=-1000;

  d_frame=0;
  int i=0;
  // set the really standard names: ENER, ENERES, VOLPRT, RLAM and FREN  
  for(; i<num_ener; i++){
    ostringstream os;
    os << "ENER[" << i+1 << "]";
    d_map.insert(MP(os.str(), i));
  }
  for(; i<num_eneres; i++){
    ostringstream os;
    os << "ENERES[" << i-num_ener+1 << "]";
    d_map.insert(MP(os.str(), i));
  }
  for(; i<num_volprt; i++){
    ostringstream os;
    os << "VOLPRT[" << i-num_eneres+1 << "]";
    d_map.insert(MP(os.str(), i));
  }
  for(; i<num_rlam; i++){
    d_map.insert(MP("RLAM", i));
  }
  for(; i<num_fren; i++){
    ostringstream os;
    os << "FREN[" << i-num_rlam+1 << "]";
    d_map.insert(MP(os.str(), i));
  }
}

int EnergyTraj::index(string s)
{
  map<string, int>::const_iterator iter=d_map.find(s);
  if(iter!=d_map.end()) return iter->second;
  
  // otherwise it might be one of the charge groups
  string sub=s.substr(0,6);
  int g_num = atoi(s.substr(s.find("[")+1, s.find("]")).c_str())-1;

  if(sub=="ENERLJ") return num_fren + 4*g_num;
  if(sub=="ENERCL") return num_fren + 4*g_num + 1;
  if(sub=="ENERRF") return num_fren + 4*g_num + 2;
  if(sub=="ENERRC") return num_fren + 4*g_num + 3;

  return unknownvariable;
}

    
double EnergyTraj::operator[](int i)
{
  if(i!=unknownvariable){
    
    if(i>=0 && i<num_fren+d_negr)
      return d_data[i];
    else if (i<0){
      
      //calculate
      int ind= -i-1;
      if(d_recalc[ind]){
	
	vector<double> v;
	for(unsigned int k=0; k< d_dep[ind].size(); k++)
	  v.push_back((*this)[d_dep[ind][k]]);
	d_e[ind].setValues(v);
	d_recalc[ind]=false;
	d_calc[ind]=d_e[ind].value();
      }
      return d_calc[ind];
    }
  }
  
  throw gromos::Exception("EnergyTraj", 
			  "Trying to access an unkwown element");
  return 0.0;
}

double EnergyTraj::operator[](string s)
{
  return (*this)[this->index(s)];
}


void EnergyTraj::addKnown(string s, string e)
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
    ostringstream os;
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
    
    gmath::Expression e(os.str());
    
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

void EnergyTraj::read_free_energy_frame(gio::Ginstream& is)
{
  if(!d_frame)
    d_data.resize(num_fren, 0.0);
  string sdum;
  while(sdum!="FREEENERGYLAMBDA") is >> sdum;
  
  // first read in the first nine elements of the energy block
  int i=0;
  for(; i< num_enerf; i++)
    is >> d_data[i];
  // rlam
  i=num_volprt;
  for(; i< num_rlam; i++)
    is >> d_data[i];
  
  // the fren block
  for(; i< num_fren; i++)
    is >> d_data[i];
  is >> sdum >> sdum;
  
  //set all things that need to be calculated to be uncalculated
  for(unsigned int i=0; i<d_recalc.size(); i++) d_recalc[i]=true;
}

  
void EnergyTraj::read_energy_frame(gio::Ginstream& is)
{
  if(!d_frame)
    d_data.resize(num_fren,0.0);
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
    d_data.resize(num_fren+4*size);
  for(int j=0; j<size; j++){
    d_data[num_fren+j]=energrp[j];
  }
  d_frame++;
  is >> sdum >> sdum;

  //set all things that need to be calculated to be uncalculated
  for(unsigned int i=0; i<d_recalc.size(); i++) d_recalc[i]=true;
  
}

string EnergyTraj::back_index(int i)
{
  map<string, int>::const_iterator iter=d_map.begin(), to=d_map.end();
  for(; iter!=to; ++iter)
    if(iter->second == i) return iter->first;
  if(i>=num_fren){
    ostringstream os;
    
    int j=i-num_fren;
    switch(j%4){
      case 0: os << "ENERLJ[" << int(j/4)+1 << "]";
	break;
      case 1: os << "ENERCL[" << int(j/4)+1 << "]";
	break;
      case 2: os << "ENERRF[" << int(j/4)+1 << "]";
	break;
      case 3: os << "ENERRC[" << int(j/4)+1 << "]";
	break;
    }
    return os.str();
  }
  
  return "unknown";
}

void EnergyTraj::write_map(ostream& os)
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

void EnergyTraj::Tokenize(const string& str,
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


  

