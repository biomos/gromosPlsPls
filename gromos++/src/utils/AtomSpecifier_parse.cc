#include <cassert>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <map>
#include <stdexcept>
#include <fstream>

#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"
#include "../gio/Ginstream.h"
#include "../bound/Boundary.h"

#include "AtomSpecifier.h"
#include "parse.h"
#include "ExpressionParser.h"

using namespace gcore;
using namespace std;
using utils::AtomSpecifier;

void AtomSpecifier::parse(std::string s, int x)
{
  std::string::size_type it = find_par(s, ';');

  if (it == std::string::npos){
    parse_single(s, x);
  }
  else{
    // std::cerr << "AS::parse\tx = " << x << std::endl;
    // std::cerr << "\tparsing " << s.substr(0, it) << std::endl;
    parse_single(s.substr(0, it), x);
    // std::cerr << "\tparsing " << s.substr(it+1, std::string::npos) << std::endl;
    parse(s.substr(it+1, std::string::npos), x);
  }
}

void AtomSpecifier::parse_single(std::string s, int x)
{
  if (s.substr(0,3) == "va("){
    std::string::size_type it = find_matching_bracket(s, '(', 3);
    parse_va(s.substr(3, it - 4), x);
  }
  else if (s.substr(0,5) == "file("){
    std::string::size_type it = find_matching_bracket(s, '(', 5);
    parse_atominfo(s.substr(5, it - 6));
  }
  else{
    std::string::size_type it = s.find(':');
    if (it == std::string::npos)
      throw(Exception("no : in AtomSpecifier"));
    
    std::vector<int> mol;
    parse_molecule(s.substr(0, it), mol, x);
    for(unsigned int i=0; i<mol.size(); ++i)
      parse_atom(mol[i], s.substr(it+1, std::string::npos), x);
  }
}

void AtomSpecifier::parse_molecule(std::string s, std::vector<int> & mol, int x)
{
  if (s=="a"){
    for(int i=0; i < d_sys->numMolecules(); ++i)
      mol.push_back(i+1);
  }
  else if(s=="s"){
    mol.push_back(0);
  }
  else if (s=="x"){
    mol.push_back(x);
  }
  else{
    // parse_mol_range(s, mol);
    // use function from parse
    parse_range(s, mol, x);
  }
}

void AtomSpecifier::parse_atom(int mol, std::string s, int x)
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

    parse_res(mol, resn, atom, x);
  }
  else{
    if (mol > 0)
      parse_atom_range(mol, 0, d_sys->mol(mol-1).numAtoms(), s, x);
    else
      parse_atom_range(mol, 0, d_sys->sol(0).numAtoms(), s, x);
  }
}

void AtomSpecifier::parse_va(std::string s, int x)
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


void AtomSpecifier::parse_atominfo(std::string s)
{
  // std::cerr << "trying to parse atominfo file" << std::endl;
  
  std::ifstream aif(s.c_str());
  if (!aif.is_open()){
    throw Exception("could not open atominfo file");
  }
  
  gio::Ginstream ai(aif);
  
  std::vector<std::string> buffer;
  ai.getblock(buffer);
  if (!buffer.size() || buffer[0] != "ATOMS"){
    std::ostringstream os;
    os << "no ATOMS block found in " << s << "!";
    throw Exception(os.str());
  }
  
  for(unsigned int i=1; i<buffer.size()-1; ++i){
    std::string s = buffer[i];
    std::string::size_type beg = s.find_first_not_of(" ");
    std::string::size_type c = s.find(':');

    int a,m;
    std::string mol_string = s.substr(beg, c-beg);
    std::istringstream is(mol_string);

    if (mol_string == "s"){
      m = 0;
    }
    else{
      if (!(is >> m)){
	std::ostringstream os;
	os << "Could not parse line: " << buffer[i] << ": trying to read mol from " << s.substr(0,c);
	throw Exception(os.str());
      }
    }
    
    is.clear();
    is.str(s.substr(c+1, std::string::npos));
    if (!(is >> a)){
      std::ostringstream os;
      os << "Could not parse line: " << buffer[i] << ": trying to read atom from " << s.substr(c+1, std::string::npos);
      throw Exception(os.str());
    }
    
    // std::cerr << "trying to add " << m << ":" << a << std::endl;
    addAtom(m-1, a-1);
  }
}


void AtomSpecifier::parse_atom_range(int mol, int beg, int end, std::string s, int x)
{
  std::map<std::string, int> var;
  var["x"] = x;
  
  ExpressionParser<int> ep;
  
  if (find_par(s, ':') != std::string::npos)
      throw Exception("Unexpected ':' token in atom set/range parsing.");

  std::string::size_type it = find_par(s, ',');
  
  if (it == std::string::npos){
    
    std::string::size_type r_it = find_par(s, '-');

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
	std::string::size_type bra = s.find_first_not_of(" ");
	int i;
	if (s[bra] == '('){
	  i = ep.parse_expression(s, var);
	  addAtom(mol-1, beg+i-1);
	}
	else{
	  // single number (or type)
	  std::istringstream is(s);
	  if(!(is >> i)){
	    if(mol > 0)
	      addType(mol-1, s, beg, end);
	    else
	      addType(mol-1, s);
	  }
	  else{
	    if ((beg + i) > end && mol > 0)
	      throw Exception("Atom out of range");	  
	    addAtom(mol-1, beg+i-1);
	  }
	}
      }
    }
    else{
      int beg, end;
      beg = ep.parse_expression(s.substr(0, r_it), var);
      end = ep.parse_expression(s.substr(r_it+1, std::string::npos), var);
      for(int i=beg; i<=end; ++i){
	addAtom(mol-1, i-1);
      }
    }
  }
  else{
    parse_atom_range(mol, beg, end, s.substr(0, it), x);
    parse_atom_range(mol, beg, end, s.substr(it+1, std::string::npos), x);
  }

}

void AtomSpecifier::parse_res(int mol, std::string res, std::string atom, int x)
{
  if (mol<0) throw Exception("No residues in solvent");

  std::map<std::string, int> var;
  var["x"] = x;
  
  ExpressionParser<int> ep;

  // std::string::size_type it = res.find(',');
  std::string::size_type it = find_par(res, ',');

  if (it == std::string::npos){
    
    // std::string::size_type r_it = res.find('-');
    std::string::size_type r_it = find_par(res, '-');

    if (r_it == std::string::npos){
      std::string::size_type bra = res.find_first_not_of(" ");
      int i;
      if (res[bra] == '('){
	i = ep.parse_expression(res, var);
	int beg, end;
	res_range(mol, i, beg, end);
	parse_atom_range(mol, beg, end, atom, x);
      }
      else{
	// single number (or type)
	std::istringstream is(res);
	int i;
	if(!(is >> i)){
	  parse_res_type(mol, res, atom, x);
	}
	else{
	  int beg, end;
	  res_range(mol, i, beg, end);
	  parse_atom_range(mol, beg, end, atom, x);
	}
      }
    }
    else{
      int beg, end;
      beg = ep.parse_expression(res.substr(0, r_it), var);
      end = ep.parse_expression
	(res.substr(r_it+1, std::string::npos),	var);

      for(int i=beg; i<=end; ++i){
	int rbeg, rend;
	res_range(mol, i, rbeg, rend);
	parse_atom_range(mol, rbeg, rend, atom, x);
      }
    }
  }
  else{
    parse_res(mol, res.substr(0, it), atom, x);
    parse_res(mol, res.substr(it+1, std::string::npos), atom, x);
  }
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

void AtomSpecifier::parse_res_type(int mol, std::string res, std::string s, int x)
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
	parse_atom_range(mol, beg, i, s, x);
	beg = i;
	resn = d_sys->mol(mol-1).topology().resNum(i);
      }
    }
    else{
      if (match){
	match = false;
	parse_atom_range(mol, beg, i, s, x);
      }
    }
  }
  if (match){
    match = false;
    parse_atom_range(mol, beg, d_sys->mol(mol-1).numAtoms(), s, x);
  }
}


