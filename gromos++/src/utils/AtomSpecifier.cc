#include <cassert>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include "AtomSpecifier.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"

using namespace gcore;
using namespace std;
using utils::AtomSpecifier;


void AtomSpecifier::parse(string s)
{
  int mol;
  int count=0;
  std::string::size_type iterator, last, colon;
  std::string token;
  

  iterator = s.find(':');
  colon=iterator;
  
  if (iterator == std::string::npos){
    if(sscanf((s.substr(0,s.length())).c_str(), "%d", &mol)!=1 &&
       s[0]!='a' && s[0]!='s')
      throw AtomSpecifier::Exception(" bad format of molecule.\n");
    if(s[0]=='a'){
      for(int m=0; m<d_sys->numMolecules(); m++){
	for(int i=0; i<d_sys->mol(m).topology().numAtoms(); i++){
	  _appendAtom(m,i);
	}
      }
    }
    else if(s[0]=='s'){
      for(int a=0; a<d_sys->sol(0).topology().numAtoms(); a++)
	addSolventType(d_sys->sol(0).topology().atom(a).name());
    }
    else{
      --mol;
      assert(mol>=0);
      for(int i=0; i<d_sys->mol(mol).topology().numAtoms(); i++){
	_appendAtom(mol,i);
      }
    }
  }
  else{
    if (sscanf((s.substr(0,iterator)).c_str(), "%d", &mol) !=1 &&
        s.substr(0,2)!="a:"&&s.substr(0,2)!="s:")
      throw AtomSpecifier::Exception(" bad format of first molecule.\n");
    int molb=0;
    int mole=d_sys->numMolecules();
    if(s.substr(0,2)!="a:"){
      molb=mol-1;
      mole=mol;
    }
    if(s.substr(0,2)=="s:"){
      molb=-1;
      mole=0;
    }
    
    for(mol=molb;mol<mole;mol++){
      
      assert(mol>=-1);
      iterator=colon;
      
      for(;true;count++){
	last=iterator;
	iterator = s.find(',',iterator+1);
	if (iterator==std::string::npos){
	  _parseAtomsHelper(s.substr(last+1, s.length()), mol);
	  break;
	} 
	_parseAtomsHelper(s.substr(last+1, iterator-last-1), mol);
      }
    }
  }
  
}
void AtomSpecifier::_parseAtomsHelper(std::string substring, int &mol)
{
  std::string::size_type iterator;
  int atomb, atome;
  if((iterator=substring.find(':')) != std::string::npos){
    if(sscanf(substring.substr(0,iterator).c_str(), "%d", &mol) != 1)
      throw AtomSpecifier::Exception(" substring: bad molecule format.\n");
    --mol;
    substring=substring.substr(iterator+1, substring.length());
  }
  if ((iterator = substring.find('-')) != std::string::npos){
    // add a range...
    if ((sscanf(substring.substr(0,iterator).c_str(), "%d", &atomb) != 1)
	|| (sscanf(substring.substr(iterator+1, substring.length()).c_str(),
		   "%d", &atome) != 1))
      throw AtomSpecifier::Exception(" substring: bad range format.\n");
    // correct for 0 indexing into arrays
    --atomb;
    --atome;
    // sanity check
    assert(atomb >= 0 && atomb < atome);
    // check whether there are enough atoms in molecule
    if (mol>=0 && d_sys->mol(mol).numAtoms() <= atome)
      throw AtomSpecifier::Exception(" not enough atoms in molecule.\n");
    
    // add the range
    for(int i=atomb; i<=atome; i++)
      _appendAtom(mol,i);
    
    return;
  }
  // adding single atom
  if (sscanf(substring.c_str(), "%d", &atomb) != 1)
    // assume that it is a name, and add the type
    addType(mol, substring);
  else{
    // correct for 0 indexing into arrays
    --atomb;
    // sanity check
    assert(atomb >= 0);
    if (mol>=0 && d_sys->mol(mol).numAtoms() <= atomb)
      throw AtomSpecifier::Exception(" not enough atoms in molecule.\n");
    _appendAtom(mol, atomb);
  }
}

void AtomSpecifier::_appendAtom(int m, int a)
{
  // check whether it is already in
  if(findAtom(m, a) == -1){
    d_atom.push_back(a);
    d_mol.push_back(m);
  }
}

bool AtomSpecifier::_checkName(int m, int a, std::string s)
{
  std::string::size_type iterator=s.find('?');
  std::string name_in_topo;
  // take in the following three lines to allow for solvent as well
  if(m<0) 
    name_in_topo=d_sys->sol(0).topology().atom(a).name().substr(0, iterator);
  else
    name_in_topo=d_sys->mol(m).topology().atom(a).name().substr(0, iterator);
  
  if (s.substr(0, iterator) == name_in_topo)
    return true;
  else 
    return false;
}

int AtomSpecifier::_expandSolvent()
{
  int nsa = d_sys->sol(0).topology().numAtoms();
  d_nsm = d_sys->sol(0).numPos() / nsa;
  
  // first remove all atoms that are in the list due to an earlier 
  // expansion. These have d_mol[i] == -2
  for(unsigned int i=0; i<d_mol.size(); i++)
    if(d_mol[i]==-2) {
      removeAtom(i); 
      i--;
    }
  // now add the atoms for every molecule
  for(int i=0; i< d_nsm; i++){
    for(unsigned int j=0; j< d_solventType.size(); j++){
      _appendAtom(-2,i*nsa+j);
    }
  }

  return d_mol.size();
}

bool AtomSpecifier::_expand()
{
  return (d_nsm != 
	  d_sys->sol(0).numPos() / d_sys->sol(0).topology().numAtoms());
}
   
AtomSpecifier::AtomSpecifier(gcore::System &sys)
{
  d_sys=&sys;
  // set d_nsm to something weird so that it will be expanded on first use
  d_nsm=-1;
}
AtomSpecifier::AtomSpecifier(gcore::System &sys, string s)
{
  d_sys=&sys;
  // set d_nsm to something weird so that it will be expanded on first use
  d_nsm=-1;
  parse(s);
}
void AtomSpecifier::setSystem(gcore::System &sys)
{
  d_sys=&sys;
  // set the number of solvent molecules to something weird so that
  // we are sure the will be correctly expanded for the new system
  d_nsm=-1;
}

int AtomSpecifier::addSpecifier(string s)
{
  parse(s);
  return d_mol.size();
}

int AtomSpecifier::addAtom(int m, int a)
{
  if(m>=d_sys->numMolecules())
    throw AtomSpecifier::Exception(
    " molecule number out of range.\n");  
  if(m>=0)
    if(a>=d_sys->mol(m).topology().numAtoms())
      throw AtomSpecifier::Exception(
      "atom number out of range.\n");
  _appendAtom(m,a);

  return d_mol.size();
}

int AtomSpecifier::addAtomStrict(int m, int a)
{
  if(m>=d_sys->numMolecules())
    throw AtomSpecifier::Exception(
    " molecule number out of range.\n");  
  if(m>=0)
    if(a>=d_sys->mol(m).topology().numAtoms())
      throw AtomSpecifier::Exception(
      "atom number out of range.\n");
   d_atom.push_back(a);
   d_mol.push_back(m);

  return d_mol.size();
}


int AtomSpecifier::addType(int m, std::string s)
{
  //loop over all atoms
  if(m<0)
    addSolventType(s);
  else{
    for(int j=0; j<d_sys->mol(m).numAtoms(); j++)
      if(_checkName(m, j, s))
	_appendAtom(m,j);
  }
  return d_mol.size();
}

int AtomSpecifier::addSolventType(std::string s)
{
  for(int j=0; j<d_sys->sol(0).topology().numAtoms(); j++){
    if(_checkName(-1,j,s)){
      int found=0;
      for(unsigned int i=0; i<d_solventType.size(); i++)
	if(d_solventType[i]==j) found =1;
      if(!found)
	d_solventType.push_back(j);
    }
  }
  _expandSolvent();
  return d_solventType.size();
}

int AtomSpecifier::addType(std::string s)
{
  //loop over all solute atoms
  for(int m=0;m<d_sys->numMolecules();m++)
    addType(m,s);
  // and do solvent
  addType(-1,s);

  return d_mol.size();
}

bool AtomSpecifier::_compare(int i, int m, int a)
{
  if(d_mol[i]==m) return d_atom[i] > a;
  if(d_mol[i]>=0 && m>=0) return d_mol[i] > m;
  if(d_mol[i]<0 && m<0) return d_atom[i] > a;
  if(d_mol[i]<0 && m>=0) return true;
  if(d_mol[i]>=0 && m<0) return false;
  return false;
}

void AtomSpecifier::sort()
{
  for(unsigned int i=1; i<d_mol.size(); i++){
    int t_m=d_mol[i], t_a=d_atom[i];
    int j=i-1;
    while ((j>=0) && _compare(j,t_m,t_a)){
      d_mol[j+1]=d_mol[j];
      d_atom[j+1]=d_atom[j];
      j--;
    }
    
    d_mol[j+1]=t_m;
    d_atom[j+1]=t_a;
  }
}

  
int AtomSpecifier::removeAtom(int m, int a)
{
  int i=findAtom(m,a);
  return removeAtom(i);
}

int AtomSpecifier::removeAtom(int i)
{
  if(i<int(d_mol.size()) && i>=0){
    vector<int>::iterator itm=d_mol.begin()+i;
    vector<int>::iterator ita=d_atom.begin()+i;
    d_mol.erase(itm);
    d_atom.erase(ita);
  }
  return d_mol.size();
}


int AtomSpecifier::findAtom(int m, int a)
{
  vector<int>::iterator itm=d_mol.begin();
  vector<int>::iterator ita=d_atom.begin();
  int counter=0;
  // a bit of a nuisance that m=-1 and m=-2 could both mean the same in this
  // function. So split up the cases
  if(m<0){
    for(; itm!=d_mol.end(); itm++, ita++, counter++)
      if(*itm<0 && *ita==a)
	return counter;
  }
  
  for(;itm!=d_mol.end(); itm++, ita++, counter++)
    if(*itm==m&&*ita==a)
      return counter;

  return -1;
}

AtomSpecifier::AtomSpecifier &AtomSpecifier::operator=(const AtomSpecifier &as)
{
  if(this != &as){
    d_mol.erase(d_mol.begin(), d_mol.end());
    d_atom.erase(d_atom.begin(), d_atom.end());
    d_solventType.erase(d_solventType.begin(), d_solventType.end());
    d_sys=as.d_sys;
    d_nsm=as.d_nsm;
    for(unsigned int i=0;i<as.d_mol.size(); i++){
      d_mol.insert(d_mol.end(), as.d_mol[i]);
      d_atom.insert(d_atom.end(), as.d_atom[i]);
    }
    for(unsigned int i=0;i<as.d_solventType.size(); i++)
      d_solventType.push_back(as.d_solventType[i]);
  }
  return *this;
  
}
AtomSpecifier::AtomSpecifier AtomSpecifier::operator+(const AtomSpecifier &as)
{
  AtomSpecifier temp(*as.d_sys);
  temp = *this;
  for(unsigned int i=0;i<as.d_mol.size();i++)
    temp.addAtom(as.d_mol[i], as.d_atom[i]);
  for(unsigned int i=0;i<as.d_solventType.size(); i++)
    temp.addSolventType(
      d_sys->sol(0).topology().atom(as.d_solventType[i]).name());
  
  return temp;
}
gmath::Vec *AtomSpecifier::coord(int i)
{
  if(_expand()) _expandSolvent();
  if(d_mol[i]<0){
    // for the solvent this test has to be done here, because we do
    // not know in advance how many water molecules we will have.
    if(d_atom[i] >= d_sys->sol(0).numPos())
      throw(AtomSpecifier::Exception("Not enough solvent in the system"));
    return &(d_sys->sol(0).pos(d_atom[i]));
  }
  else
    return &(d_sys->mol(d_mol[i]).pos(d_atom[i]));
}
  
std::string AtomSpecifier::name(int i)
{
  if(_expand()) _expandSolvent();
  if(d_mol[i] < 0){
    int num=(d_atom[i]%d_sys->sol(0).topology().numAtoms());
    return d_sys->sol(0).topology().atom(num).name();
  }
  else
    return d_sys->mol(d_mol[i]).topology().atom(d_atom[i]).name();
} 
int AtomSpecifier::iac(int i)
{
  if(_expand()) _expandSolvent();
  if(d_mol[i] < 0){
    int num=(d_atom[i]%d_sys->sol(0).topology().numAtoms());
    return d_sys->sol(0).topology().atom(num).iac();
  }
  else
    return d_sys->mol(d_mol[i]).topology().atom(d_atom[i]).iac();
} 
double AtomSpecifier::charge(int i)
{
  if(_expand()) _expandSolvent();
  if(d_mol[i] < 0){
    int num=(d_atom[i]%d_sys->sol(0).topology().numAtoms());
    return d_sys->sol(0).topology().atom(num).charge();
  }
  else
    return d_sys->mol(d_mol[i]).topology().atom(d_atom[i]).charge();
}

int AtomSpecifier::mol(int i)
{
  if(_expand()) _expandSolvent();
  if(i>=int(d_mol.size()))
    throw AtomSpecifier::Exception("Trying to access non-existing element");
  return d_mol[i];
}
int AtomSpecifier::atom(int i)
{
  if(_expand()) _expandSolvent();
  if(i>=int(d_atom.size())){
    throw AtomSpecifier::Exception("Trying to access non-existing element");
  }
  return d_atom[i];
}

int AtomSpecifier::size()
{
  if(_expand()) _expandSolvent();
  return d_mol.size();
}

void AtomSpecifier::clear()
{
  d_mol.resize(0);
  d_atom.resize(0);
  d_solventType.resize(0);
  d_nsm=-1;
}
std::vector<std::string> AtomSpecifier::toString()
{
 
    
  std::vector<std::string> s;
  ostringstream os;
  if(size()){
    sort();
    if(d_mol[0]<0) os << "s"; else os << d_mol[0] + 1;
    os << ":" << d_atom[0] + 1;
    
    int lastatom=d_atom[0];
    string::size_type fPos=0, lPos=0;
  
    for(unsigned int i=1; i< d_mol.size(); i++){
      if(d_mol[i]==d_mol[i-1]){
	if(d_atom[i] != d_atom[i-1] + 1){
	  if(d_atom[i-1] > lastatom)
	    os << "-" << d_atom[i-1] + 1;
	  lPos= os.str().length();
	  s.push_back(os.str().substr(fPos, lPos-fPos));
	  fPos=lPos;
	  if(d_mol[i] < 0) os << "s"; else os << d_mol[i] + 1;
	  os << ":" << d_atom[i]+1;
	  lastatom = d_atom[i];
	}
      }    
      else{
	if(d_atom[i-1] > lastatom)
	  os << "-" << d_atom[i-1] + 1;
	lPos= os.str().length();
	s.push_back(os.str().substr(fPos, lPos));
	fPos=lPos;
	if(d_mol[i]<0) os << "s";
	else os << d_mol[i] + 1;
	os << ":" << d_atom[i] + 1;
	lastatom = d_atom[i];
      }
    }
    // has the last atom been closed?
    if(d_atom[d_atom.size()-1] != lastatom){
      os << "-" << d_atom[d_atom.size()-1]+1;
    }
    lPos= os.str().length();
    s.push_back(os.str().substr(fPos, lPos));
  }
  
  return s;
}
