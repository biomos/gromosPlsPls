#include "InNoe.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../rest/VirtualAtom.h"
#include "../rest/Neighbours.h"
#include <strstream>
#include <vector>
#include <algorithm>
#include <cmath>


using utils::InNoe_i;
using utils::InNoe;
using namespace gcore;
using namespace rest;


class InNoe_i{
  friend InNoe;
  
  const System &d_sys;
  vector<VirtualAtom*> d_at[2];
  vector<double> d_dist;
  InNoe_i(const System &sys): d_sys(sys), d_at(), d_dist(){}
  ~InNoe_i(){}
};



InNoe::InNoe(const System  &sys, const string &line):
  d_this(new InNoe_i(sys))
{
  char buff[80];
  string sat[2];
  double d;
  
  // get the two atoms
  istrstream is(line.c_str());
  is>>buff;
  sat[0]=buff;
  is>>buff;
  sat[1]=buff;

  // get all distances
  if(!is.good())
    throw Exception("InNoe specified by \n" + line + "\ncontains no values for the InNoe!");
  while(is>>d){
    d_this->d_dist.push_back(d);
  }

  // parse the atoms to InNoes
  for(int k=0;k<2;++k){

    // What type of InNoe is it?
    /* It can be none or 'E' for an explicit InNoe, 
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
      while(at > (atNum+=sys.mol(mol).numAtoms())) ++mol;
      at-=atNum+sys.mol(mol).numAtoms();
      

      switch(type){

      case 'E':
	// explicit restraint
	if(i!=sat[k].size()-1)
	  throw Exception("Inconsistent input line:\n"+line);
	d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 5));
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
	  while(at > (atNum+=sys.mol(mol).numAtoms())) ++mol;
	  at-=atNum+sys.mol(mol).numAtoms();
	  
	  // push second one...
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 5));

	  break;

	case 2:
	  // stereospecific CH2, type 4 l and r
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 4));
	  d_this->d_at[k].push_back(new VirtualAtom(sys,mol,at, 4, 'l'));

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
}

  

const double InNoe::corr_dist(const int i)const
{
  double cd=dist(i);
  for(int k=0;k<2;k++){
    switch(InNoe_[0].type(k)){ // all InNoe[i] HAVE to be of the same type...
    case 3:
      cd+=.09;
      break;
    case 5:
      cd*=pow(3.0,1.0/3.0);
      cd+=.03;
      break;
    case 6:
      cd+=.22;
      break;
    case 7:
      cd+=.21;
      break;
    }
  }
  
  return cd;
  
}


const string InNoe::info()const
{
  ostrstream os;
  os << "# InNoe consisting of "<< InNoe_.size() << " possible distances:\n";
  for (unsigned int i=0;i<InNoe_.size(); i++)
    os << "# " << ' ' << InNoe_[i].info() << endl;
  os << '\0';
  
  return os.str();
  
}

