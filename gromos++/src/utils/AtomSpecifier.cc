#include <iostream>
#include <string>
#include "AtomSpecifier.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"

using namespace gcore;
using namespace std;
using utils::AtomSpecifier;

static void parse(AtomSpecifier &as, string s)
{
  string m;
  int mol=-1, ab=-1, ae=-1;
  
  int i, max=s.length();
  for(i=0;s[i]!=':'&&i!=max;i++)
    m+=s[i];
  mol=atoi(m.c_str())-1;
  m="";
  if(mol<0)
    throw AtomSpecifier::Exception(
    " invalid atom-specifier.\nSyntax: <mol>[:<first atom>[-<last atom>]]\n");
  if(mol>=as.sys()->numMolecules())
    throw AtomSpecifier::Exception(
    " molecule number out of range.\n");  
  if(i==max){
    ab=0;
    ae=as.sys()->mol(mol).topology().numAtoms();
  }
  else {
    for(i++;s[i]!='-'&&i!=max;i++)
      m+=s[i];
    ab=atoi(m.c_str())-1;
    m="";
    if(ab<0)
      throw AtomSpecifier::Exception(
      " invalid atom-specifier.\nSyntax: <mol>[:<first atom>[-<last atom>]]\n");
    if(ab>=as.sys()->mol(mol).topology().numAtoms())
      throw AtomSpecifier::Exception(
      "first atom number out of range (>numAtoms)\n");
    if(i==max) {
      ae=ab+1;
    }
    else{
      for(i++;i!=max;i++)
        m+=s[i];
      ae=atoi(m.c_str());
      if(ae<0)
        throw AtomSpecifier::Exception(
     " invalid atom-specifier.\nSyntax: <mol>[:<first atom>[-<last atom>]]\n");
      if(ae>as.sys()->mol(mol).topology().numAtoms()||ae<=ab)
        throw AtomSpecifier::Exception(
         "last atom number out of range (< first atom or > numAtoms)\n");
    }
    
  }
  // now we fill up the vectors of as
  for(int j=ab;j<ae;j++){
    as.mol()->insert(as.mol()->end(), mol);
    as.atom()->insert(as.atom()->end(), j);
  }
  
}

AtomSpecifier::AtomSpecifier(gcore::System &sys)
{
  d_sys=&sys;
}
AtomSpecifier::AtomSpecifier(gcore::System &sys, string s)
{
  d_sys=&sys;
  parse(*this, s);
}
int AtomSpecifier::addSpecifier(string s)
{
  parse(*this, s);
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


