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

void AtomSpecifier::parse(std::string s)
{
  std::string::size_type it = s.find(';');
  
  if (it == std::string::npos){
    parse_single(s);
  }
  else{
    parse_single(s.substr(0, it));
    parse(s.substr(it+1, std::string::npos));
  }
}

void AtomSpecifier::parse_single(std::string s)
{
  if (s.substr(0,3) == "va("){
    std::string::size_type it = find_matching_bracket(s, '(', 3);
    parse_va(s.substr(3, it - 4));
  }
  else{
    std::string::size_type it = s.find(':');
    if (it == std::string::npos)
      throw(Exception("no : in AtomSpecifier"));
    
    std::vector<int> mol;
    parse_molecule(s.substr(0, it), mol);
    for(unsigned int i=0; i<mol.size(); ++i)
      parse_atom(mol[i], s.substr(it+1, std::string::npos));
  }
}

void AtomSpecifier::parse_molecule(std::string s, std::vector<int> & mol)
{
  if (s=="a"){
    for(int i=0; i < d_sys->numMolecules(); ++i)
      mol.push_back(i+1);
  }
  else if(s=="s"){
    mol.push_back(0);
  }
  else{
    parse_mol_range(s, mol);
  }
}

void AtomSpecifier::parse_atom(int mol, std::string s)
{
  if (s.substr(0, 4) == "res("){
    std::string::size_type ket = find_matching_bracket(s, '(', 4);
    if (ket == std::string::npos)
      throw Exception("Residue: end bracket missing");

    std::string res = s.substr(4, ket - 5);
    
    std::string::size_type sep = res.find(':');
    if (sep == std::string::npos)
      throw Exception("No atoms for residue given");
  
    std::string resn = res.substr(0, sep);
    std::string atom = res.substr(sep+1, std::string::npos);

    parse_res(mol, resn, atom);
  }
  else{
    if (mol > 0)
      parse_atom_range(mol, 0, d_sys->mol(mol-1).numAtoms(), s);
    else
      parse_atom_range(mol, 0, d_sys->sol(0).topology().numAtoms(), s);
  }
}

void AtomSpecifier::parse_va(std::string s)
{
  std::string::size_type it = s.find(',');
  if (it == std::string::npos)
    throw Exception(" Virtual Atom: wrong format. Should be va(type, as)");
  
  std::string t = s.substr(0, it);
  VirtualAtom::virtual_type vt;

  // try the types
  if (t == "com") vt = VirtualAtom::COM;
  else if (t == "cog") vt = VirtualAtom::COG;
  else{
    std::istringstream is(t);
    int i;
    if (!(is >> i))
      throw Exception(" Virtual Atom: type not recognised: " + t);
    vt = VirtualAtom::virtual_type(i);
  }
  d_specatom.push_back
    (new VirtualSpecAtom(*d_sys, 
			 s.substr(it+1, std::string::npos),
			 vt));
}

void AtomSpecifier::parse_mol_range(std::string s, std::vector<int> & range)
{
  std::string::size_type it = s.find(',');
  if (it == std::string::npos){
    
    std::string::size_type r_it = s.find('-');
    if (r_it == std::string::npos){
      // single number (or type)
      std::istringstream is(s);
      int i;
      if(!(is >> i))
	throw Exception("type in molecule not allowed");
      else{
	range.push_back(i);
      }
    }
    else{
      int beg, end;
      std::istringstream is(s.substr(0, r_it));
      if (!(is >> beg))
	throw Exception("range: begin is not a number");
      is.clear();
      is.str(s.substr(r_it+1, std::string::npos));
      if (!(is >> end))
	throw Exception("range: end is not a number");
      for(int i=beg; i<=end; ++i)
	range.push_back(i);
    }
  }
  else{
    parse_mol_range(s.substr(0, it), range);
    parse_mol_range(s.substr(it+1, std::string::npos), range);
  }
}

void AtomSpecifier::parse_atom_range(int mol, int beg, int end, std::string s)
{
  std::string::size_type it = s.find(',');
  if (it == std::string::npos){
    
    std::string::size_type r_it = s.find('-');
    if (r_it == std::string::npos){

      if (s == "a"){
	if(mol >0)
	  for(int i=beg; i<end; ++i)
	    addAtom(mol-1, i);
	else
	  for(int i=0; i<d_sys->sol(0).topology().numAtoms(); ++i)
	    addSolventType(d_sys->sol(0).topology().atom(i).name());
      }
      else{
	// single number (or type)
	std::istringstream is(s);
	int i;
	if(!(is >> i)){
	  if(mol > 0)
	    addType(mol-1, s, beg, end);
	  else
	    addType(mol-1, s);
	}
	else{
	  if ((beg + i) > end)
	    throw Exception("Atom out of range");	  
	  addAtom(mol-1, beg+i-1);
	}
      }
    }
    else{
      int beg, end;
      std::istringstream is(s.substr(0, r_it));
      if (!(is >> beg))
	throw Exception("range: begin is not a number");
      is.clear();
      is.str(s.substr(r_it+1, std::string::npos));
      if (!(is >> end))
	throw Exception("range: end is not a number");
      for(int i=beg; i<=end; ++i){
	addAtom(mol-1, i-1);
      }
    }
  }
  else{
    parse_atom_range(mol, beg, end, s.substr(0, it));
    parse_atom_range(mol, beg, end, s.substr(it+1, std::string::npos));
  }

}

void AtomSpecifier::parse_res(int mol, std::string res, std::string atom)
{
  if (mol<0) throw Exception("No residues in solvent");

  std::string::size_type it = res.find(',');
  if (it == std::string::npos){
    
    std::string::size_type r_it = res.find('-');
    if (r_it == std::string::npos){
      // single number (or type)
      std::istringstream is(res);
      int i;
      if(!(is >> i)){
	parse_res_type(mol, res, atom);
      }
      else{
	int beg, end;
	res_range(mol, i, beg, end);
	parse_atom_range(mol, beg, end, atom);
      }
    }
    else{
      int beg, end;
      std::istringstream is(res.substr(0, r_it));
      if (!(is >> beg))
	throw Exception("range: begin is not a number");
      is.clear();
      is.str(res.substr(r_it+1, std::string::npos));
      if (!(is >> end))
	throw Exception("range: end is not a number");
      for(int i=beg; i<=end; ++i){
	int rbeg, rend;
	res_range(mol, i, rbeg, rend);
	parse_atom_range(mol, rbeg, rend, atom);
      }
    }
  }
  else{
    parse_res(mol, res.substr(0, it), atom);
    parse_res(mol, res.substr(it+1, std::string::npos), atom);
  }
}

std::string::size_type AtomSpecifier::find_matching_bracket
(
 std::string s,
 char bra,
 std::string::size_type it)
{
  char ket;
  if (bra == '(') ket = ')';
  else if (bra == '[') ket = ']';
  else if (bra == '{') ket = '}';
  else
    throw Exception("Bracket not recognised");
  
  int level = 1;
  for( ; it < s.length() && level != 0; ++it){
    if (s[it] == bra) ++level;
    else if (s[it] == ket) --level;
  }
  
  if (level) return std::string::npos;
  else return it;
}

void AtomSpecifier::res_range(int mol, int res, int &beg, int &end)
{
  beg = d_sys->mol(mol-1).numAtoms();
  end = 0;
  --res;
  
  for(int i=0; i<d_sys->mol(mol-1).numAtoms(); ++i){
    if (d_sys->mol(mol-1).topology().resNum(i) == res){
      if (i < beg) beg = i;
      if (i > end) end = i;
    }
  }
  ++end;
}

void AtomSpecifier::parse_res_type(int mol, std::string res, std::string s)
{
  int beg = 0;
  bool match=false;
  int resn = 0;

  for(int i=0; i<d_sys->mol(mol-1).numAtoms(); ++i){
    
    if (_checkResName(mol-1, i, res)){
      if(!match){
	beg = i;
	match = true;
	resn = d_sys->mol(mol-1).topology().resNum(i);
      }
      else if(resn != d_sys->mol(mol-1).topology().resNum(i)){
	parse_atom_range(mol, beg, i, s);
	beg = i;
	resn = d_sys->mol(mol-1).topology().resNum(i);
      }
    }
    else{
      if (match){
	match = false;
	parse_atom_range(mol, beg, i, s);
      }
    }
  }
  if (match){
    match = false;
    parse_atom_range(mol, beg, d_sys->mol(mol-1).numAtoms(), s);
  }
}
