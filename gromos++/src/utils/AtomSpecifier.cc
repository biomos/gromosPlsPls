#include <cassert>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>

#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"

#include "AtomSpecifier.h"

using namespace gcore;
using namespace std;
using utils::AtomSpecifier;

std::string utils::SpecAtom::toString()const
{
  std::ostringstream os;
  os << d_mol << ":" << d_atom;
  return os.str();
}

void AtomSpecifier::parse(std::string s)
{
  // take care of virtual atoms
  std::string::size_type v_it = s.find('(');
  if (v_it != std::string::npos){
    // hooray
    if ((v_it == 2 && s.substr(v_it - 2, 3) != "va(") ||
	(v_it > 2 && s.substr(v_it - 3, 4) != ",va("))
      throw AtomSpecifier::Exception(" wrong format for virtual atom");

    int level = 1;
    std::string::size_type v_end = v_it + 1;
    while(level && v_end < s.length()){
      if (s[v_end] == '(') ++level;
      if (s[v_end] == ')') --level;
      ++v_end;
    }
    if (level || (v_end < s.length() && s[v_end] != ','))
      throw AtomSpecifier::Exception(" missing bracket or comma in virtual atom");
    
    // let's split up
    if (v_it > 2){
      std::string s1 = s.substr(0, v_it - 3);
      parse(s1);
    }
    
    std::string s2 = s.substr(v_it + 1, v_end - v_it - 2);
    _parseVirtualAtom(s2);

    if (v_end + 1 < s.length()){
      std::string s3 = s.substr(v_end + 1, std::string::npos);
      parse(s3);
    }
    
    return;
  }

  std::string::size_type it = s.find(':');
  std::string ss;
  
  int mol, molb = 0, mole = 0;
  
  if (it == std::string::npos)
    _parseWholeMolecule(s);

  while(s != ""){

    it = s.find(',');
    if (it != std::string::npos){
      ss = s.substr(0, it);
      s = s.substr(it + 1, std::string::npos);
    }
    else{
      ss = s;
      s = "";
    }

    // check for a molecule (":")
    std::string::size_type col = ss.find(':');
    if (col != std::string::npos){
      std::string mols = ss.substr(0, col);
      ss = ss.substr(col + 1, std::string::npos);

      if (mols == "a"){
	molb = 0;
	mole = d_sys->numMolecules();
      }
      else if (mols == "s"){
	molb = -1;
	mole = 0;
      }
      else{
	std::istringstream imol(mols);
	if (!(imol >> mol))
	  throw AtomSpecifier::Exception(" bad format of molecule: "+s+"\n");
	molb = mol - 1;
	mole = mol;
      }
    }

    for(mol = molb; mol < mole; ++mol){
      assert(mol>=-1);
      _parseAtomsHelper(ss, mol);
    }

  } // all "," seperated parts done

}

void AtomSpecifier::_parseWholeMolecule(string s)
{
  // we're adding a whole molecule...
  std::istringstream iss(s);
  int mol;
  std::string trail;
  
  if((iss >> mol)==0 && s!="a" && s!="s") 
    throw AtomSpecifier::Exception(" bad format of molecule: "+s+"\n");
  if((iss >> trail)!=0)
    throw AtomSpecifier::Exception(" trailing data in molecule: "+s+
				   "\ndon't understand " +trail+ "\n");
    
  if(s=="a"){
    // easy: that's ALL atoms (of all molecules)
    for(int m=0; m<d_sys->numMolecules(); m++){
      addMolecule(mol);
    }
  }

  else if(s=="s"){
    // all solvent atoms (of solvent 0...)
    for(int a=0; a<d_sys->sol(0).topology().numAtoms(); a++)
      addSolventType(d_sys->sol(0).topology().atom(a).name());
  }
  else{
    --mol;
    if(mol<0)
      throw AtomSpecifier::Exception(" molecule numbers should be > 0\n");
    if(mol>=d_sys->numMolecules())
      throw AtomSpecifier::Exception(" not enough molecules in system: " 
				     +s+"\n");

    addMolecule(mol);

  }
}

void AtomSpecifier::_parseAtomsHelper(std::string substring, int &mol)
{
  std::string::size_type iterator;
  int atomb, atome;

  if((iterator=substring.find(':')) != std::string::npos){
    throw AtomSpecifier::Exception(" substring: bad atom format (contains :).\n");
  }

  if ((iterator = substring.find('-')) != std::string::npos){
    // add a range...
    std::istringstream ir(substring.substr(0, iterator));
    if (!(ir >> atomb))
      throw AtomSpecifier::Exception(" substring: bad range format (begin)\n"+substring);
    ir.clear();
    ir.str(substring.substr(iterator + 1, std::string::npos));
    if (!(ir >> atome))
      throw AtomSpecifier::Exception(" substring: bad range format (end)\n"+substring);

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
  std::istringstream ir(substring.substr(0, iterator));
  if (!(ir >> atomb)){
    // assume that it is a name, and add the type
    addType(mol, substring);
  }
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

void AtomSpecifier::_parseVirtualAtom(std::string s)
{
  std::string::size_type it = s.find(',');
  if (it == std::string::npos)
    throw AtomSpecifier::Exception(" Virtual Atom: wrong format. Should be va(type, as)");
  
  std::string t = s.substr(0, it);
  VirtualAtom::virtual_type vt;

  // try the types
  if (t == "com") vt = VirtualAtom::COM;
  else if (t == "cog") vt = VirtualAtom::COG;
  else{
    std::istringstream is(t);
    int i;
    if (!(is >> i))
      throw AtomSpecifier::Exception(" Virtual Atom: type not recognised: " + t);
    vt = VirtualAtom::virtual_type(i);
  }
  
  d_specatom.push_back(new VirtualSpecAtom(*d_sys, s.substr(it+1, std::string::npos), vt));

}

void AtomSpecifier::_appendAtom(int m, int a)
{
  // check whether it is already in
  if(findAtom(m, a) == -1){
    if (m < 0)
      d_specatom.push_back(new SolventSpecAtom(*d_sys, m, a));
    else
      d_specatom.push_back(new SpecAtom(*d_sys, m, a));
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
  for(unsigned int i=0; i<d_specatom.size(); ++i){
    if(d_specatom[i]->mol()==-2) {
      removeAtom(i); 
      i--;
    }
  }
  
  // now add the atoms for every molecule
  for(int i=0; i < d_nsm; ++i){
    for(unsigned int j=0; j< d_solventType.size(); ++j){
      _appendAtom(-2,i*nsa+d_solventType[j]);
    }
  }

  return d_specatom.size();
}

bool AtomSpecifier::_expand()
{
  return (d_nsm != 
	  d_sys->sol(0).numPos() / d_sys->sol(0).topology().numAtoms());
}
   
AtomSpecifier::AtomSpecifier(gcore::System &sys)
{
  d_sys = &sys;
  // set d_nsm to something weird so that it will be expanded on first use
  d_nsm = -1;
}

AtomSpecifier::AtomSpecifier(gcore::System &sys, string s)
{
  d_sys = &sys;
  // set d_nsm to something weird so that it will be expanded on first use
  d_nsm = -1;
  parse(s);
}

AtomSpecifier::~AtomSpecifier()
{
  // delete all atoms from the spec
  // (this is slightly undefined behaviour, as the deleted pointers
  //  are left in the vector (for clear))
  for(std::vector<SpecAtom *>::iterator it = d_specatom.begin(),
	to = d_specatom.end(); it != to; ++it){
    delete *it;
  }
}

void AtomSpecifier::setSystem(gcore::System &sys)
{
  d_sys = &sys;
  // set the number of solvent molecules to something weird so that
  // we are sure the will be correctly expanded for the new system
  d_nsm = -1;

  // need to change the systems in the SpecAtom vector
  for(unsigned int i=0; i<d_specatom.size(); ++i)
    d_specatom[i]->setSystem(sys);
}

int AtomSpecifier::addSpecifier(string s)
{
  parse(s);
  return d_specatom.size();
}

int AtomSpecifier::addAtom(int m, int a)
{
  if(m >= d_sys->numMolecules())
    throw AtomSpecifier::Exception(" molecule number out of range.\n");  
  if(m >= 0)
    if(a >= d_sys->mol(m).topology().numAtoms())
      throw AtomSpecifier::Exception(" atom number out of range.\n");
  _appendAtom(m,a);
  
  return d_specatom.size();
}

int AtomSpecifier::addGromosAtom(int a)
{
  int m=0;
  
  while(a >= d_sys->mol(m).numAtoms()){
    a -= d_sys->mol(m).numAtoms();
    m++;
    if(m >= d_sys->numMolecules()){
      m=-1;
      break;
    }
  }

  if(m >= 0)
    if(a >= d_sys->mol(m).topology().numAtoms())
      throw AtomSpecifier::Exception(" atom number out of range.\n");

  _appendAtom(m,a);
  return d_specatom.size();
}

int AtomSpecifier::addMolecule(int m)
{
  if(m >= d_sys->numMolecules())
    throw AtomSpecifier::Exception(" molecule number out of range.\n");  
  
  for(int i=0; i < d_sys->mol(m).numAtoms(); ++i)
    _appendAtom(m, i);

  return d_specatom.size();

}

int AtomSpecifier::addAtomStrict(int m, int a)
{
  if(m >= d_sys->numMolecules())
    throw AtomSpecifier::Exception(" molecule number out of range.\n");  
  if(m >= 0){
    if(a >= d_sys->mol(m).topology().numAtoms())
      throw AtomSpecifier::Exception(" atom number out of range.\n");
    
    d_specatom.push_back(new SpecAtom(*d_sys, a, m));
    
  }
  else{
    d_specatom.push_back(new SolventSpecAtom(*d_sys, a, m));
  }

  return d_specatom.size();
}


int AtomSpecifier::addType(int m, std::string s)
{
  //loop over all atoms
  if(m<0)
    addSolventType(s);
  else{
    for(int j=0; j<d_sys->mol(m).numAtoms(); ++j)
      if(_checkName(m, j, s))
	_appendAtom(m,j);
  }
  return d_specatom.size();
}

int AtomSpecifier::addSolventType(std::string s)
{
  for(int j=0; j < d_sys->sol(0).topology().numAtoms(); ++j){
    if(_checkName(-1,j,s)){
      int found=0;
      for(unsigned int i = 0; i < d_solventType.size(); ++i)
	if(d_solventType[i]==j) found = 1;
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
  for(int m = 0; m < d_sys->numMolecules(); ++m)
    addType(m, s);
  // and do solvent
  addType(-1, s);

  return d_specatom.size();
}

bool AtomSpecifier::_compare(int i, int m, int a)
{
  if(d_specatom[i]->mol() == m) return d_specatom[i]->atom() > a;
  if(d_specatom[i]->mol() >= 0 && m >= 0) return d_specatom[i]->mol() > m;
  if(d_specatom[i]->mol() <  0 && m <  0) return d_specatom[i]->atom() > a;
  if(d_specatom[i]->mol() <  0 && m >= 0) return true;
  if(d_specatom[i]->mol() >= 0 && m <  0) return false;
  return false;
}

void AtomSpecifier::sort()
{
  for(unsigned int i=1; i < d_specatom.size(); ++i){
    SpecAtom * t = d_specatom[i];
    int j=i-1;
    while ((j >= 0) && _compare(j, t->mol(), t->atom())){
      d_specatom[j+1] = d_specatom[j];
      j--;
    }
    d_specatom[j+1] = t;
  }
}
  
int AtomSpecifier::removeAtom(int m, int a)
{
  int i=findAtom(m,a);
  return removeAtom(i);
}

int AtomSpecifier::removeAtom(int i)
{
  if(i < int(d_specatom.size()) && i >= 0){
    vector<SpecAtom *>::iterator it=d_specatom.begin() + i;
    SpecAtom * d = *it;
    
    d_specatom.erase(it);
    delete d;
  }
  return d_specatom.size();
}

int AtomSpecifier::findAtom(int m, int a)
{
  vector<SpecAtom *>::iterator it=d_specatom.begin(),
    it_to = d_specatom.end();
  
  int counter=0;
  // a bit of a nuisance that m=-1 and m=-2 could both mean the same in this
  // function. So split up the cases
  if(m<0){
    for( ; it != it_to; ++it, ++counter)
      if((*it)->mol() < 0 && (*it)->atom() == a)
	return counter;
  }
  
  for( ; it != it_to; ++it, ++counter)
    if((*it)->mol() == m && (*it)->atom() == a)
      return counter;

  return -1;
}

/**
 * copy constructor.
 */
AtomSpecifier::AtomSpecifier(utils::AtomSpecifier const & as)
{
  if (this != &as){
    clear();
    
    d_sys=as.d_sys;
    d_nsm=as.d_nsm;

    std::vector<SpecAtom *>::const_iterator
      it = as.d_specatom.begin(),
      to = as.d_specatom.end();
    for( ; it != to; ++it){
      d_specatom.push_back((*it)->clone());
    }

    d_solventType = as.d_solventType;
  }
}

AtomSpecifier::AtomSpecifier &AtomSpecifier::operator=(const AtomSpecifier &as)
{
  if(this != &as){

    clear();
    
    d_sys=as.d_sys;
    d_nsm=as.d_nsm;

    for(unsigned int i=0; i < as.d_specatom.size(); ++i){
      d_specatom.push_back(as.atom()[i]->clone());
    }

    d_solventType = as.d_solventType;
  }

  return *this;
}

AtomSpecifier::AtomSpecifier AtomSpecifier::operator+(const AtomSpecifier &as)
{
  AtomSpecifier temp(*as.d_sys);
  temp = *this;

  for(unsigned int i = 0; i < as.d_specatom.size(); ++i)
    temp.addAtom(as.d_specatom[i]->mol(), as.d_specatom[i]->atom());
  
  for(unsigned int i=0;i<as.d_solventType.size(); i++)
    temp.addSolventType(d_sys->sol(0).topology().atom(as.d_solventType[i]).name());
  
  // if copy construction works here, why not use it in the first statement???
  return temp;
}

gmath::Vec *AtomSpecifier::coord(int i)
{
  if(_expand()) _expandSolvent();
  return &d_specatom[i]->pos();
}

gmath::Vec & AtomSpecifier::pos(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->pos();
}
  
std::string AtomSpecifier::name(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->name();
} 

int AtomSpecifier::iac(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->iac();
} 

double AtomSpecifier::charge(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->charge();
}
double AtomSpecifier::mass(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->mass();
}

std::string AtomSpecifier::resname(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->resname();
}

int AtomSpecifier::resnum(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->resnum();
}

int AtomSpecifier::mol(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->mol();
}

int AtomSpecifier::atom(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->atom();
}

int AtomSpecifier::size()
{
  if(_expand()) _expandSolvent();
  return d_specatom.size();
}

void AtomSpecifier::clear()
{
  std::vector<SpecAtom *>::iterator
    it = d_specatom.begin(),
    to = d_specatom.end();
  
  for( ; it != to; ++it)
    delete *it;
  
  d_specatom.clear();
  d_solventType.resize(0);
  d_nsm=-1;
}

std::vector<std::string> AtomSpecifier::toString()
{

  std::vector<std::string> s;
  ostringstream os;
  int m = -999, a_last = -999;
  bool in_range = false, first = true;
  
  for(unsigned int i = 0; i < d_specatom.size(); ++i){
    
    // virtual atom
    if (d_specatom[i]->type() == spec_virtual){
      if (in_range){
	in_range = false;
	os << a_last + 1;
      }
      os << "," << d_specatom[i]->toString();
      m = -999;
      a_last = -999;
      continue;
    }
    
    // new molecule?
    if (d_specatom[i]->mol() != m){
      m = d_specatom[i]->mol();
      
      if (in_range){
	in_range = false;
	os << "-" << a_last + 1;
      }

      a_last = d_specatom[i]->atom();

      if (!first){
	os << ",";
	if (m < 0) os << "s";
	else os << m + 1;
	
	os << ":" << a_last + 1;
      }
      else{
	if (m < 0) os << "s";
	else os << m+1;
	
	os << ":" << a_last + 1;
	first = false;
      }
    }
    else if (a_last == d_specatom[i]->atom()-1){
      a_last = d_specatom[i]->atom();
      in_range = true;
    }
    else if (in_range){
      in_range = false;
      os << "-" << a_last + 1 << ",";
      a_last = d_specatom[i]->atom();
      os << a_last + 1;
    }
    else{
      a_last = d_specatom[i]->atom();
      os << "," << a_last + 1;
    }
  }

  // but the why???
  s.push_back(os.str());
  return s;
  
  /*
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
  */
  
}
