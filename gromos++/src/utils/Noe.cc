#include <cassert>
//#include <strstream>

#include "Noe.h"
#include "../gio/StringTokenizer.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "VirtualAtom.h"
#include "Neighbours.h"
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace gio;
using namespace std;
using namespace gcore;
using utils::Noe_i;
using utils::Noe;


class Noe_i{
  friend class Noe;

  const System &d_sys;
  std::vector<VirtualAtom*> d_at[2];
  std::vector<double> d_dist;
  std::vector<double> cor;
  std::vector<int> cortype;
  int d_num;  
  Noe_i(const System &sys): d_sys(sys), d_at(), d_dist(){}
  ~Noe_i(){}
};



Noe::Noe(const System  &sys, const string &line):d_this(new Noe_i(sys))
{
  //set default correction values
  d_this->cor.push_back(0.0);
  d_this->cor.push_back(0.0);
  d_this->cor.push_back(0.0);
  d_this->cor.push_back(0.0);
  
  d_this->cortype.push_back(3);
  d_this->cortype.push_back(5);
  d_this->cortype.push_back(6);
  d_this->cortype.push_back(7);



  //parse line
  char buff[80];
  string sat[2];
  double d;
  
  // get the two atoms
  istringstream is(line.c_str());
  is>>buff;
  sat[0]=buff;
  is>>buff;
  sat[1]=buff;

  // get all distances
  if(!is.good())
    throw Exception("Noe specified by \n" + line + "\ncontains no values for the Noe!");
  while(is>>d){
    d_this->d_dist.push_back(d);
  }
  sort(d_this->d_dist.begin(), d_this->d_dist.end());
  
  // parse the atoms to Noes
  for(int k=0;k<2;++k){

    // What type of Noe is it?
    /* It can be none or 'E' for an explicit Noe, 
       'N' for a non-stereospecific CHn group,
       'S' for a stereospecific CHn group,
       'A' for an aromatic H or a rotating aromatic ring */
    
    string::size_type i=sat[k].find_first_of("ENSA");
    
    // check for valid and unique specification
    if(i<sat[k].size()){ 
      string::size_type j=sat[k].find_first_not_of("0123456789");
      if(j < sat[k].size() && j!=i)
	throw Exception("Inconsistent input line:\n"+line);
      j=sat[k].find_last_not_of("0123456789");
      if(j < sat[k].size() && j!=i)
	throw Exception("Inconsistent input line:\n"+line);
    }
    else{ // default for 'E'
      i=sat[k].size();
      sat[k]+='E';
    }
    
    
    // Now do the parsing
    try{

      char type=sat[k][i];
      int at=atoi(sat[k].substr(0,i).c_str())-1;
      int mol=0, atNum=0;

      // parse into mol and atom rather than high atom nr.
      while(at >= (atNum+=sys.mol(mol).numAtoms())){
	++mol;
	if(mol >= sys.numMolecules())
	  throw Exception("Atom number too high in input line:\n"+line);
      }
      at-=atNum-sys.mol(mol).numAtoms();
      

      switch(type){

      case 'E':
	// explicit restraint
	if(i!=sat[k].size()-1)
	  throw Exception("Inconsistent input line:\n"+line);
	d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 0));
	break;
	
	
      case 'N': 
	// non-stereospecific CHn

	// is N the last letter? -> this is only possibility
	if(i!=sat[k].size()-1)
	  throw Exception("Inconsistent input line:\n"+line);

	switch(Neighbours(sys,mol,at).size()){
	case 1:
	  // CH3, type 5
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 5));
	  break;
	case 2:
	  // CH2, type 3
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 3));
	  break;
	case 3:
	  // non-stereospecific rotating methyl groups (Val, Leu), type 6
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 6));
	  break;
	default:
	  throw Exception("Inconsistent input line:\n"+line);
	}

	break;
	
      case 'S':
	switch(Neighbours(sys,mol,at).size()){
	case 1:
	  // stereospecific rotating methyls (Val, Leu), 2 x type 5

	  // is there a second one?
	  if(i==sat[k].size()-1)
	    throw Exception("Inconsistent input line:\n"+line);

	  // push first one...
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 5));
	  
	  // parse second one...
	  at=atoi(sat[k].substr(i+1).c_str())-1;
	  mol=0;
	  atNum=0;
	  while(at >= (atNum+=sys.mol(mol).numAtoms())) ++mol;
	  at-=atNum-sys.mol(mol).numAtoms();
	  
	  // push second one...
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 5));

	  break;

	case 2:
	  // stereospecific CH2, type 4 l and r
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 4));
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 4, 1));

	  break;
	  
	case 3:
	  // stereospecific CH1, type 1
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 1));
	  break;


	default:
	  throw Exception("Inconsistent input line:\n"+line);
	}
	
	break;

      case 'A':
	
	// is A the last letter? -> this is only possibility
	if(i!=sat[k].size()-1)
	  throw Exception("Inconsistent input line:\n"+line);

	switch(Neighbours(sys,mol,at).size()){
	case 3:
	  // rotating ring, type 7
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 7));
	  break;
	  
	case 2:
	  // CH1 (aromatic, old ff), type 2
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 2));
	  break;
	  
	default:
	  throw Exception("Inconsistent input line:\n"+line);
	}
	break;
	
      default:
	throw Exception("Inconsistent input line:\n"+line);
      } /* end switch(type)*/
    } /* end try */
    catch(VirtualAtom::Exception &e){
      string mess=string(e.what())+"\nInput line: "+line;
      throw Exception(mess);
    }
  } /* end for */

  d_this->d_num =  d_this->d_at[0].size()*d_this->d_at[1].size();

  if( int(d_this->d_dist.size()) > d_this->d_num)
    throw Exception("Too many distances specified in input line:\n" +line);
  
} // Noe::Noe(const System  &sys, const string &line):d_this(new Noe_i(sys))

//implementation of the new NOE program/constructor, taking a PROADR input
Noe::Noe(const System  &sys, const string &line, double dish, double disc):d_this(new Noe_i(sys))
{

  //d_this->Dish = dish;
  //d_this->Disc = disc;
  //parse line into Tokens...
  StringTokenizer stok(line, " ");
  std::vector<std::string> tokens = stok.tokenize();
  if (tokens.size() < 11) throw Exception("At least 11 input digits are expected! Check input!\n");

  //smack in the distance
    d_this->d_dist.push_back(atof(tokens[10].c_str()));

    //smack in the atomic noe information
  int offset = 0;
  for(int k=0;k<2;++k) {

    int at = atoi(tokens[0+offset].c_str())-1;
    int type = atoi(tokens[4+offset].c_str());
    int mol=0, atNum=0;
    
    // parse into mol and atom rather than high atom nr.
    while(at >= (atNum+=sys.mol(mol).numAtoms())){
      ++mol;
      if(mol >= sys.numMolecules())
	throw Exception("Atom number too high in input line:\n"+line);
    }

    // atoms are always bound
    atNum-=sys.mol(mol).numAtoms();
    
    at-=atNum;
    
    int config[4];
    
    config[0] = atoi(tokens[0+offset].c_str())-1-atNum;
    config[1] = atoi(tokens[1+offset].c_str())-1-atNum;
    config[2] = atoi(tokens[2+offset].c_str())-1-atNum;
    config[3] = atoi(tokens[3+offset].c_str())-1-atNum;
    
    
    d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, type, config, dish, disc));
    //arse
    
    offset = 5;
  }
  d_this->d_num =  d_this->d_at[0].size()*d_this->d_at[1].size();
  
}

double Noe::distance(int i)const{
  assert(i<d_this->d_num);
  int at[2];

  at[1] = i % d_this->d_at[1].size();
  at[0] = (i-at[1])/d_this->d_at[1].size();

  return ( d_this->d_at[0][at[0]]->pos() -
	   d_this->d_at[1][at[1]]->pos() ).abs();

}  

string Noe::info(int i)const{
  assert(i<d_this->d_num);
  ostringstream os;

  int at[2];

  at[1] = i % d_this->d_at[1].size();
  at[0] = (i-at[1])/d_this->d_at[1].size();

  for (int j=0 ; j < 2; ++j){
    int atom = d_this->d_at[j][at[j]]->operator[](0);
    int mol = d_this->d_at[j][at[j]]->mol();
    int type = d_this->d_at[j][at[j]]->type();

    int resNum = d_this->d_sys.mol(mol).topology().resNum(atom);
    string resName = d_this->d_sys.mol(mol).topology().resName(resNum);
    string atName = d_this->d_sys.mol(mol).topology().atom(atom).name();

     
    //std::ostrstream oss;
     ostringstream oss;
    if (type == 4) {
     oss << " Type 4 Noe! Atoms: "
         << ((d_this->d_at[j][at[j]]->operator[](0))+1) << " " 
         << ((d_this->d_at[j][at[j]]->operator[](1))+1) << " " 
         << ((d_this->d_at[j][at[j]]->operator[](2))+1) << " " 
         << ((d_this->d_at[j][at[j]]->operator[](3))+1);// << ends;
    }

    std::string typefour(oss.str());

    //os.freeze(false);  // Unfreeze so buffer is deleted.

    //cout << typefour << endl;
  
    
    switch(type){
    case 3:
    case 5:
    case 8:
    case 6:
      atName[0] = 'Q';
      break;
    case 4:
      atName += typefour;
    case 1:
    case 2:
      atName[0] = 'H';
      break;
    }
      
      

    os.setf(ios::right, ios::adjustfield);
    os << setw(5) << resNum+1;
    os.setf(ios::left, ios::adjustfield);
    os <<  setw(5) << resName.c_str() << setw(5) << atName.c_str();
  }

  return string(os.str());
}


int Noe::numDistances()const{
  return d_this->d_num;
}
int Noe::numReferences()const{
  return d_this->d_dist.size();
}


string Noe::distRes(int i)const{

  assert(i<d_this->d_num);
  int at[2];

  at[1] = i % d_this->d_at[1].size();
  at[0] = (i-at[1])/d_this->d_at[1].size();


  ostringstream ss;
  ss.setf(ios::right, ios::adjustfield);
  ss.setf(ios::fixed, ios::floatfield);
  ss.precision(3);

  for(int j=0; j< 2; ++j){
   for (int k=0;k<4;++k){
      int att = d_this->d_at[j][at[j]]->operator[](k);
      int mol = d_this->d_at[j][at[j]]->mol();
      int offset = 1;
      for(int j=0;j<mol;++j)
	offset += d_this->d_sys.mol(j).numAtoms();
      if(att == -1)
	ss << setw(5) << 0;
      else
	ss << setw(5) << att + offset;
    }
    ss << setw(3) << d_this->d_at[j][at[j]]->type() % 7;
  }
  

  // do some ugly, ugly crap
   if (i < int(d_this->d_dist.size())){
  ss << setw(10) << correctedReference(i); 
  }
 else if (i == int(d_this->d_dist.size())){
   ss << setw(10) << correctedReference(i-1);}
  else if (i == 2){
   ss << setw(10) << correctedReference(i-2);}
  else if (i == 3){
   ss << setw(10) << correctedReference(i-3);}
  
  ss << setw(10) << 1.0;
 
  return string(ss.str());
}

double Noe::reference(int i)const{
  assert( i < int(d_this->d_dist.size()));
  return d_this->d_dist[i];
}

double Noe::correctedReference(int i)const{
  assert( i < int(d_this->d_dist.size()));
  double cd=d_this->d_dist[i];

  //now come several multiplicity corrections
  //see: Neuhaus D. et.al. J.Bio.NMR, 8 (1996) 292-310

  // for type 5, the experimental distance has to be multiplied by
  // 3^(1/6)
  for(int k=0;k<2;k++){
      if(d_this->d_at[k][0]->type()==5) {
	cd*=pow(3.0,1.0/6.0); 
      }
  // for type 3, the experimental distance has to be multiplied by
  // 2^(1/6)
      else if(d_this->d_at[k][0]->type()==3) {
        cd*=pow(2.0,1.0/6.0);
      }
  // for type 6, the experimental distance has to be multiplied by
  // 6^(1/6)
      else if(d_this->d_at[k][0]->type()==6) {
        cd*=pow(6.0,1.0/6.0);
      }
  // for type 7, the experimental distance has to be multiplied by
  // 2^(1/6)
      else if(d_this->d_at[k][0]->type()==7) {
        cd*=pow(2.0,1.0/6.0);
      }
  }

  //parse the actual read in corrections from the correction file

  for(int k=0;k<2;k++){
    switch(d_this->d_at[k][0]->type()){ 

    case 3:
      cd+=d_this->cor[0];
      break;
    case 5:
      cd+=d_this->cor[1];
      break;
    case 6:
      cd+=d_this->cor[2];
      break;
    case 7:
      cd+=d_this->cor[3];
      break;
    }
  }
  
  return cd;
  
}

double Noe::correction(int type) {

    if ((type < 3) || (type == 4) || (type > 7))
      throw Exception("GROMOS Noe type not known:" +type);
  
  double t=0;
  for (int i=0; i < int (d_this->cortype.size()); ++i){
    if (d_this->cortype[i] == type) {
      t = d_this->cor[i];}
  }
    
  return t;
}
     

void Noe::setcorrection(int type, double correction) {
  if ((type < 3) || (type == 4) || (type > 7))
      throw Exception("GROMOS Noe type not known:" +type);
  
  for (int i=0; i < int (d_this->cortype.size()); ++i){
    if (d_this->cortype[i] == type) {
      d_this->cor[i] = correction;}
  }
}





