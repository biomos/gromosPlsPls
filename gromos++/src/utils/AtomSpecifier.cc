#include <iostream>
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
       s[0]!='a')
      throw AtomSpecifier::Exception(" bad format of molecule.\n");
    if(s[0]=='a'){
      for(int m=0; m<d_sys->numMolecules(); m++){
	for(int i=0; i<d_sys->mol(m).topology().numAtoms(); i++){
	  d_atom.insert(d_atom.end(), i);
	  d_mol.insert(d_mol.end(), m);
	}
      }
    }
    else{
      --mol;
      assert(mol>=0);
      for(int i=0; i<d_sys->mol(mol).topology().numAtoms(); i++){
	d_atom.insert(d_atom.end(), i);
	d_mol.insert(d_mol.end(), mol);
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
  if ((iterator = substring.find('-')) != std::string::npos)
    {
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
          {
            d_atom.insert(d_atom.end(), i);
            d_mol.insert(d_mol.end(), mol);
          }
        return;
      }
    // adding single atom
    if (sscanf(substring.c_str(), "%d", &atomb) != 1)
      throw AtomSpecifier::Exception(" substring: bad atom format.\n");

    // correct for 0 indexing into arrays
    --atomb;
    // sanity check
    assert(atomb >= 0);
    if (mol>=0 && d_sys->mol(mol).numAtoms() <= atomb)
      throw AtomSpecifier::Exception(" not enough atoms in molecule.\n");
    
    d_atom.insert(d_atom.end(), atomb);
    d_mol.insert(d_mol.end(), mol);
       
}

AtomSpecifier::AtomSpecifier(gcore::System &sys)
{
  d_sys=&sys;
}
AtomSpecifier::AtomSpecifier(gcore::System &sys, string s)
{
  d_sys=&sys;
  parse(s);
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
  d_mol.insert(d_mol.end(), m);
  d_atom.insert(d_atom.end(), a);
  
  return d_mol.size();
}
int AtomSpecifier::addType(std::string s)
{
  //loop over all atoms
  for(int i=0;i<d_sys->numMolecules();i++)
    for(int j=0;j<d_sys->mol(i).topology().numAtoms(); j++)
      if((s==d_sys->mol(i).topology().atom(j).name())||s=="ALL")
	addAtom(i,j);
  return d_mol.size();
}

int AtomSpecifier::removeAtom(int m, int a)
{
  vector<int>::iterator itm=d_mol.begin();
  vector<int>::iterator ita=d_atom.begin();
  for(;itm!=d_mol.end();itm++){
    if(*itm==m&&*ita==a) {
      d_mol.erase(itm);
      d_atom.erase(ita);
      break;
      
    }
    ita++;
  }
  
  return d_mol.size();  
}

int AtomSpecifier::findAtom(int m, int a)
{
  vector<int>::iterator itm=d_mol.begin();
  vector<int>::iterator ita=d_atom.begin();
  int counter=0;
  
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
    
    d_sys=as.d_sys;
    for(unsigned int i=0;i<as.d_mol.size(); i++){
      d_mol.insert(d_mol.end(), as.d_mol[i]);
      d_atom.insert(d_atom.end(), as.d_atom[i]);
    }
  }
  return *this;
  
}
AtomSpecifier::AtomSpecifier AtomSpecifier::operator+(const AtomSpecifier &as)
{
  AtomSpecifier temp(*as.d_sys);
  temp = *this;
  for(unsigned int i=0;i<as.d_mol.size();i++)
    temp.addAtom(as.d_mol[i], as.d_atom[i]);
  
  return temp;
}
gmath::Vec *AtomSpecifier::coord(int i)
{
  if(d_mol[i]<0){
    // for the solvent this test has to be done here, because we do
    // not know in advance how many water molecules we will have.
    if(d_atom[i]>=d_sys->sol(0).numCoords())
      throw(AtomSpecifier::Exception("Not enough solvent in the system"));
    
    return &(d_sys->sol(0).pos(d_atom[i]));
  }
  else
    return &(d_sys->mol(d_mol[i]).pos(d_atom[i]));
}
  
std::string AtomSpecifier::name(int i)
{
   if(d_mol[i]<0){
    // for the solvent this test has to be done here, because we do
    // not know in advance how many water molecules we will have.
    if(d_atom[i]>=d_sys->sol(0).numCoords())
      throw(AtomSpecifier::Exception("Not enough solvent in the system"));
    int num=(d_atom[i]%d_sys->sol(0).topology().numAtoms());
    return d_sys->sol(0).topology().atom(num).name();
    
  }
  else
    return d_sys->mol(d_mol[i]).topology().atom(d_atom[i]).name();
} 
