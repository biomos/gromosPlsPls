// 	$Id$	
  
//---Property Class-----------------------------------

#include "PropertyContainer.h"
#include <iostream>
#include <math.h>
#include <stdio.h>
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"

using namespace gcore;
using namespace std;
using namespace utils;

namespace utils
{
   Property::Property(gcore::System &sys)
   {
     d_type = "Property";
     REQUIREDARGUMENTS=0;
     d_ubound = MAXFLOAT;
     d_lbound = -MAXFLOAT;
     d_zvalue = 0;
     
     d_sys = &sys;
   }

  Property::Property()
  {
    // only for child classes
    d_sys = NULL;
  }
  
  Property::~Property()
  {
  }

  void Property::parse(int count, std::string arguments[])
  {
    // default implementation
    // takes care of the first four arguments
    // is called directly on the child object, therefore not virtual
  
    if (count < REQUIREDARGUMENTS)
      throw Property::Exception(" too few arguments.\n");
  
    if (count >4) count = 4; // for switch()

    switch (count){
    case 4:
      {
	if (sscanf(arguments[3].c_str(),"%f",&d_ubound) != 1)
	  throw Property::Exception(" upper boundary: bad format.\n");
      }
    case 3:
      {
	if (sscanf(arguments[2].c_str(),"%f", &d_lbound) != 1)
	  throw Property::Exception(" lower boundary: bad format.\n");
      }
    case 2:
      {
	if (sscanf(arguments[1].c_str(), "%f", &d_zvalue) != 1)
	  throw Property::Exception(" zero value: bad format.\n");
      }
    case 1:
      {
	parseAtoms(arguments[0]);
	// error handling within
	break;
      }
    case 0:
      {
	// called without arguments
	throw Property::Exception(" no property type specified (no arguments at all).\n");
      }
    }
    
  }
  
  void Property::parseAtoms(std::string atoms)
  {
    int mol;
    int count = 0;
    std::string::size_type iterator, last;
    std::string token;
    
    // a first molecule is required!
    // at least in this (standard) implementation
    iterator = atoms.find(':');
    if (iterator == std::string::npos)
      throw Property::Exception(" atoms: bad format.\n");

    if (sscanf((atoms.substr(0,iterator)).c_str(),"%d", &mol) != 1)
      throw Property::Exception(" atoms: bad format of first molecule.\n");

    // arrays are zero based:
    --mol;
    assert(mol>=0);
    // tokenize the rest of the string at commas
    for(; true; count++)
      {
	last = iterator;
	iterator = atoms.find(',',iterator+1);
	if (iterator == std::string::npos)
	  {
	    // add last substring
	    _parseAtomsHelper(atoms.substr(last+1, atoms.length()), mol);
	    count++;
	    break;
	  }
	// add substring
	_parseAtomsHelper(atoms.substr(last+1, iterator-last-1), mol);
      }
  }

  void Property::_parseAtomsHelper(std::string substring, int &mol)
  {
    // substring is of the following format:
    // [mol:]atom[-atom]

    // if mol is not present, then the mol from last token will be used...
    std::string::size_type iterator;
    int atomb, atome;
    
    // check whether there is a mol:
    if ((iterator = substring.find(':')) != std::string::npos)
      {
	// a new molecule...
	if (sscanf(substring.substr(0,iterator).c_str(), "%d", &mol) != 1)
	  throw Property::Exception(" substring: bad molecule format.\n");
	
	// indexing begins with 0
	--mol;
	// update substring (remove the '<mol>:')
	substring = substring.substr(iterator+1, substring.length());
      }
    
    // check for a range
    if ((iterator = substring.find('-')) != std::string::npos)
      {
	// add a range...
	if ((sscanf(substring.substr(0,iterator).c_str(), "%d", &atomb) != 1)
	    || (sscanf(substring.substr(iterator+1, substring.length()).c_str(),
		       "%d", &atome) != 1))
	  throw Property::Exception(" substring: bad range format.\n");
	// correct for 0 indexing into arrays
	--atomb;
	--atome;
	// sanity check
	assert(atomb >= 0 && atomb < atome);
	// check whether there are enough atoms in molecule
	if (d_sys->mol(mol).numAtoms() <= atome)
	  throw Property::Exception(" not enough atoms in molecule.\n");
	
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
      throw Property::Exception(" substring: bad atom format.\n");

    // correct for 0 indexing into arrays
    --atomb;
    // sanity check
    assert(atomb >= 0);
    if (d_sys->mol(mol).numAtoms() <= atomb)
      throw Property::Exception(" not enough atoms in molecule.\n");
    
    d_atom.insert(d_atom.end(), atomb);
    d_mol.insert(d_mol.end(), mol);
       
  }

  double Property::calc()
  {
    throw Property::Exception(" calculation (for general property) not implemented.\n");

    return 0;
  }

  std::string Property::checkBounds()
  {
    char b[100];
    std::string s = "";
    if (d_value > d_ubound || d_value < d_lbound)
      {
	sprintf(b,"# Boundary Violation: %s %f\n", toTitle().c_str(), d_value);
	s = b;
      }
    return s;
  }
  

  std::string Property::toTitle()
  { 
    return d_type + ": general";
  }

  std::string Property::toString()
  {
    return d_type + ": general";
  }

  std::string Property::average()
  {
    throw Property::Exception(" averaging not implemented.\n");
  }

  std::ostream &operator<<(std::ostream &os, Property &p)
  {
    os << p.toString();
    return os;
  }
  

  //---DistanceProperty Class------------------------------------

  DistanceProperty::DistanceProperty(gcore::System &sys) :
    Property(sys)
  {
    d_type = "Distance";
    REQUIREDARGUMENTS = 1;
    
    d_average = 0;
    d_zrmsd = 0;
    d_count = 0;
   
  }
  
  DistanceProperty::~DistanceProperty()
  {
  }

  void DistanceProperty::parse(int count, std::string arguments[])
  {
    Property::parse(count, arguments);

    // for a distance, we should just have to atoms here...
    if (d_atom.size() != 2)
      throw DistanceProperty::Exception("wrong number of atoms for a distance.\n");
  }

  double DistanceProperty::calc()
  {
    gmath::Vec tmp = ((d_sys->mol(d_mol[0])).pos(d_atom[0])-
		      (d_sys->mol(d_mol[1])).pos(d_atom[1])); 

    // save the value for a boundary check
    // and for printing!
    d_value = tmp.abs();
    // averaging
    d_average += d_value;
    d_zrmsd += pow(d_value - d_zvalue, 2);

    ++d_count;
    
    return d_value;
  }
  
  std::string DistanceProperty::toString()
  {
    char b[100];
    sprintf(b, "%f", d_value);
    std::string s = b;
    return s;
  }


  std::string DistanceProperty::toTitle()
    {
      char b[100];
      std::string s;
      sprintf(b, "d%%%d:%d-%d:%d", d_mol[0]+1, d_atom[0]+1, d_mol[1]+1, d_atom[1]+1);
      s = b;
      return s;
    }
  
  std::string DistanceProperty::average()
  {
    // return <average> <rmsd from z-value>
    char b[100];
    std::string s;
    double z = d_zrmsd / d_count;
    z = sqrt(z);
    sprintf(b, "%f\t\t%f", d_average / d_count, z);
    s = b;
    return s;
  }

  //---AngleProperty Class------------------------------------------------------------

  AngleProperty::AngleProperty(gcore::System &sys) :
    Property(sys)
  {
    d_type = "Angle";
    REQUIREDARGUMENTS = 1;
    
    d_average = 0;
    d_zrmsd = 0;
    d_count = 0;
  }
  
  AngleProperty::~AngleProperty()
  {
    // nothing to do...
  }
  
  void AngleProperty::parse(int count, std::string arguments[])
  {
    Property::parse(count, arguments);
    
    // it's an angle, therefore 3 atoms
    if (d_atom.size() != 3)
      throw AngleProperty::Exception("wrong number of atoms for an angle.\n");
  }
  
  double AngleProperty::calc()
  {
    gmath::Vec tmpA = (d_sys->mol(d_mol[0]).pos(d_atom[0])
		      -d_sys->mol(d_mol[1]).pos(d_atom[1])); 
    gmath::Vec tmpB = (d_sys->mol(d_mol[2]).pos(d_atom[2])
		      -d_sys->mol(d_mol[1]).pos(d_atom[1]));
    d_value = acos((tmpA.dot(tmpB))/(tmpA.abs()*tmpB.abs()))*180/pi;
    
    d_average += d_value;
    d_zrmsd += pow(d_value - d_zvalue, 2);

    ++d_count;
    return d_value;
  }
  
  std::string AngleProperty::toString()
  {
    char b[100];
    sprintf(b, "%f", d_value);
    std::string s = b;
    return s;
  }
  
  std::string AngleProperty::toTitle()
  {
    char b[100];
    std::string s;
    sprintf(b, "a%%%d:%d-%d:%d-%d:%d", d_mol[0]+1, d_atom[0]+1, 
	    d_mol[1]+1, d_atom[1]+1, d_mol[2]+1, d_atom[2]+1);
    s = b;
    return s;
  }
  
  std::string AngleProperty::average()
  {
    char b[100];
    std::string s;
    double z = sqrt(d_zrmsd / d_count);
    sprintf(b, "%f\t\t%f", d_average / d_count, z);
    s = b;
    return s;
  }
  

  //---TorsionProperty Class------------------------------------

  TorsionProperty::TorsionProperty(gcore::System &sys) :
    Property(sys)
  {
    d_type = "Torsion";
    REQUIREDARGUMENTS = 1;
    
    d_average = 0;
    d_zrmsd = 0;
    d_count = 0;
  }
  
  TorsionProperty::~TorsionProperty()
  {
  }
  
  void TorsionProperty::parse(int count, std::string arguments[])
  {
    Property::parse(count, arguments);
    
    // it's a torsion, therefore 4 atoms needed
    if (d_atom.size() != 4)
      throw TorsionProperty::Exception(
	      " wrong number of atoms for torsion.\n");
  }
  
  double TorsionProperty::calc()
  {
    gmath::Vec tmpA = (d_sys->mol(d_mol[0]).pos(d_atom[0])
		      -d_sys->mol(d_mol[1]).pos(d_atom[1]));
    gmath::Vec tmpB = (d_sys->mol(d_mol[3]).pos(d_atom[3])
		      -d_sys->mol(d_mol[2]).pos(d_atom[2]));
    gmath::Vec tmpC = (d_sys->mol(d_mol[2]).pos(d_atom[2])
		      -d_sys->mol(d_mol[1]).pos(d_atom[1]));

    gmath::Vec p1 = tmpA.cross(tmpC);
    gmath::Vec p2 = tmpB.cross(tmpC);

    double cosphi = ((p1.dot(p2))/(p1.abs()*p2.abs()));

    d_value = acos(cosphi)*180/pi;     

    gmath::Vec p3 = p1.cross(p2);
    if (p3.dot(tmpC)<0)
      d_value = 360 - d_value;   
    
    d_average += d_value;
    d_zrmsd += pow(d_value-d_zvalue, 2);
    
    ++d_count;
    return d_value;
  }

  std::string TorsionProperty::toString()
  {
    char b[100];
    sprintf(b, "%f", d_value);
    std::string s = b;
    return s;
  }
  
  std::string TorsionProperty::toTitle()
  {
    char b[100];
    std::string s;
    sprintf(b, "t%%%d:%d-%d:%d-%d:%d-%d:%d", d_mol[0]+1, d_atom[0]+1,
	    d_mol[1]+1, d_atom[1]+1, d_mol[2]+1, d_atom[2]+1,
	    d_mol[3]+1, d_atom[3]+1);
    s = b;
    return s;
  }
  
  std::string TorsionProperty::average()
  {
    char b[100];
    std::string s;
    double z = sqrt(d_zrmsd / d_count);
    
    sprintf(b, "%f\t\t%f", d_average / d_count, z);
    s = b;
    return s;
  }
  
}















