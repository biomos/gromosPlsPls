// 	$Id$	
  
//---Property Class-----------------------------------

#include <cassert>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <stdio.h>

#include "../gmath/Vec.h"

#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Improper.h"
#include "../gcore/Dihedral.h"
#include "../utils/AtomSpecifier.h"

#include "PropertyContainer.h"

using namespace gcore;
using namespace std;
using namespace utils;

namespace utils
{
   Property::Property(gcore::System &sys)
     : d_atom(sys)
   {
     d_type = "Property";
     REQUIREDARGUMENTS=0;
     d_ubound = MAXFLOAT;
     d_lbound = -MAXFLOAT;
     d_zvalue = 0;
     
     d_sys = &sys;
   }

  /*
  Property::Property()
  {
    // only for child classes
    d_sys = NULL;
  }
  */
  
  Property::~Property()
  {
  }

  void Property::parse(int count, std::string arguments[])
  {
    // default implementation
    // takes care of the first four arguments
  
    if (count < REQUIREDARGUMENTS)
      throw Property::Exception(" too few arguments.\n");
  
    if (count >4) count = 4; // for switch()

    switch (count){
    case 4:
      {
	if (sscanf(arguments[3].c_str(),"%lf",&d_ubound) != 1)
	  throw Property::Exception(" upper boundary: bad format.\n");
      }
    case 3:
      {
	if (sscanf(arguments[2].c_str(),"%lf", &d_lbound) != 1)
	  throw Property::Exception(" lower boundary: bad format.\n");
      }
    case 2:
      {
	if (sscanf(arguments[1].c_str(), "%lf", &d_zvalue) != 1)
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
    d_atom.addSpecifier(atoms);
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
    std::ostringstream os;
    os << d_type << "%";
    bool first = true;
    
    for(int i=0; i < atoms().size(); ++i){
      if (first) first = false;
      else os << "-";
      os << atoms().mol(i)+1 << ":" << atoms().atom(i)+1;
    }
    
    return os.str();
  }

  std::string Property::toString()
  {
    std::ostringstream os;
    os << getValue();
    return os.str();
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
  
  int Property::getTopologyType(gcore::System const &sys)
  {
    // standard implementation
    // check whether all atoms belong to the same molecule
    if (atoms().size()){
      int the_mol = atoms().mol(0);

      for(int m=1; m < atoms().size(); ++m)
	if (atoms().mol(m) != the_mol) return -1;
    
      return findTopologyType(sys.mol(the_mol).topology());
    }
    else return -1;
  }

  int Property::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    return -1;
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
    if (atoms().size() != 2)
      throw DistanceProperty::Exception("wrong number of atoms for a distance.\n");
  }

  double DistanceProperty::calc()
  {
    gmath::Vec tmp = atoms().pos(0) - atoms().pos(1);
    // save the value for a boundary check
    // and for printing!
    d_value = tmp.abs();
    // averaging
    d_average += d_value;
    d_zrmsd += pow(d_value - d_zvalue, 2);

    ++d_count;
    
    return d_value;
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

  int DistanceProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    int a, b;

    if(atoms().atom(0) < atoms().atom(1)){
      a = atoms().atom(0);
      b = atoms().atom(1);
    }
    else{
      a = atoms().atom(1);
      b = atoms().atom(0);
    }
    BondIterator bi(mol_topo);
    while(bi)
      if(bi()[0]==a && bi()[1]==b) return bi().type();
      else ++bi;
    
    return -1;
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
    if (atoms().size() != 3)
      throw AngleProperty::Exception("wrong number of atoms for an angle.\n");
  }
  
  double AngleProperty::calc()
  {
    gmath::Vec tmpA = atoms().pos(0) - atoms().pos(1);
    gmath::Vec tmpB = atoms().pos(2) - atoms().pos(1);

    d_value = acos((tmpA.dot(tmpB))/(tmpA.abs()*tmpB.abs()))*180/M_PI;
    
    d_average += d_value;
    d_zrmsd += pow(d_value - d_zvalue, 2);

    ++d_count;
    return d_value;
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
  

  int AngleProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    int a, b, c;

    if(atoms().atom(0) < atoms().atom(2)){
      a = atoms().atom(0);
      b = atoms().atom(1);
      c = atoms().atom(2);
    }
    else{
      a = atoms().atom(2);
      b = atoms().atom(1);
      c = atoms().atom(0);
    }
    AngleIterator ai(mol_topo);
    while(ai)
      if(ai()[0]==a && ai()[1]==b && ai()[2] == c) return ai().type();
      else ++ai;
    
    return -1;
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
    gmath::Vec tmpA = atoms().pos(0) - atoms().pos(1);
    gmath::Vec tmpB = atoms().pos(3) - atoms().pos(2);
    gmath::Vec tmpC = atoms().pos(2) - atoms().pos(1);

    gmath::Vec p1 = tmpA.cross(tmpC);
    gmath::Vec p2 = tmpB.cross(tmpC);

    double cosphi = ((p1.dot(p2))/(p1.abs()*p2.abs()));

    if(cosphi > 1.0) cosphi=1.0;
    if(cosphi <-1.0) cosphi=-1.0;
    d_value = acos(cosphi)*180/M_PI;     

    gmath::Vec p3 = p1.cross(p2);
    if (p3.dot(tmpC)<0)
      d_value = 360 - d_value;   
    
    d_average += d_value;
    d_zrmsd += pow(d_value-d_zvalue, 2);
    
    ++d_count;
    return d_value;
  }

  std::string TorsionProperty::average()
  {
    std::ostringstream os;
    double z = sqrt(d_zrmsd / d_count);
    
    os << std::setw(20) << d_average / d_count
       << std::setw(20) << z;

    return os.str();
  }
  
  int TorsionProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    int a, b, c, d;

    if(atoms().atom(1) < atoms().atom(2)){
      a = atoms().atom(0);
      b = atoms().atom(1);
      c = atoms().atom(2);
      d = atoms().atom(3);
    }
    else{
      a = atoms().atom(3);
      b = atoms().atom(2);
      c = atoms().atom(1);
      d = atoms().atom(0);
    }

    // assuming it is a dihedral...
    DihedralIterator di(mol_topo);
    while(di)
      if(di()[0]==a && di()[1]==b && di()[2]==c && di()[3]==d) 
	return di().type();
      else ++di;
    
    return -1;
  }


  //---OrderProperty Class------------------------------------------------------------

  OrderProperty::OrderProperty(gcore::System &sys) :
    Property(sys),
    d_average(0),
    d_zrmsd(0),
    d_count(0),
    d_axis(0)
  {
    // the first one we read ourselves...
    REQUIREDARGUMENTS = 1;
    d_type = "Order";
    
  }
  
  OrderProperty::~OrderProperty()
  {
    // nothing to do...
  }
  
  void OrderProperty::parse(int count, std::string arguments[])
  {
    if (arguments[0] == "x") d_axis = Vec(1,0,0);
    else if (arguments[0] == "y") d_axis = Vec(0,1,0);
    else if (arguments[0] == "z") d_axis = Vec(0,0,1);
    else
      throw OrderProperty::Exception("axis specification wrong.\n");

    Property::parse(count - 1, &arguments[1]);
    
    // it's an angle, therefore 3 atoms
    if (d_atom.size() != 2)
      throw OrderProperty::Exception("wrong number of atoms for an order property.\n");
  }
  
  double OrderProperty::calc()
  {
    gmath::Vec tmpA = atoms().pos(0) - atoms().pos(1);

    d_value = acos((tmpA.dot(d_axis))/(tmpA.abs()*d_axis.abs()))*180/M_PI;
    
    d_average += d_value;
    d_zrmsd += pow(d_value - d_zvalue, 2);

    ++d_count;
    return d_value;
  }
  
  std::string OrderProperty::toTitle()
  {
    std::ostringstream ss;
    ss << "o%";
    if (d_axis[0]) ss << "x";
    else if (d_axis[1]) ss << "y";
    else ss << "z";
    
    ss << "%" << atoms().mol(0)+1 << ":" << atoms().atom(0)+1
       << "-" << atoms().mol(1)+1 << ":" << atoms().atom(1)+1;
    
    return ss.str();
  }
  
  std::string OrderProperty::average()
  {
    std::ostringstream os;
    double z = sqrt(d_zrmsd / d_count);

    os << std::setw(20) << d_average / d_count
       << std::setw(20) <<  z;
    return os.str();
  }
  

  int OrderProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    // not implemented
    
    return -1;
  }

  //---OrderParamProperty Class------------------------------------------------------------

  OrderParamProperty::OrderParamProperty(gcore::System &sys) :
    Property(sys),
    d_average(0),
    d_zrmsd(0),
    d_count(0),
    d_axis(0)
  {
    // the first one we read ourselves...
    REQUIREDARGUMENTS = 1;
    d_type = "OrderParam";
    
  }
  
  OrderParamProperty::~OrderParamProperty()
  {
    // nothing to do...
  }
  
  void OrderParamProperty::parse(int count, std::string arguments[])
  {
    if (arguments[0] == "x") d_axis = Vec(1,0,0);
    else if (arguments[0] == "y") d_axis = Vec(0,1,0);
    else if (arguments[0] == "z") d_axis = Vec(0,0,1);
    else
      throw OrderParamProperty::Exception("axis specification wrong.\n");

    Property::parse(count - 1, &arguments[1]);
    
    // it's an angle, therefore 3 atoms
    if (d_atom.size() != 2)
      throw OrderParamProperty::Exception("wrong number of atoms for an order property.\n");
  }
  
  double OrderParamProperty::calc()
  {
    gmath::Vec tmpA = atoms().pos(0) - atoms().pos(1);

    const double cosa = tmpA.dot(d_axis)/(tmpA.abs()*d_axis.abs());
    d_value = 0.5 * (3 * cosa * cosa - 1);
    
    d_average += d_value;
    d_zrmsd += pow(d_value - d_zvalue, 2);

    ++d_count;
    return d_value;
  }
  
  std::string OrderParamProperty::toTitle()
  {

    std::ostringstream ss;
    ss << "op%";
    if (d_axis[0]) ss << "x";
    else if (d_axis[1]) ss << "y";
    else ss << "z";
    
    ss << "%" << atoms().mol(0)+1 << ":" << atoms().atom(0)+1
       << "-" << atoms().mol(1)+1 << ":" << atoms().atom(1)+1;
    
    return ss.str();
  }
  
  std::string OrderParamProperty::average()
  {
    std::ostringstream os;

    double z = sqrt(d_zrmsd / d_count);

    os << std::setw(20) << d_average / d_count
       << std::setw(20) << z;

    return os.str();
  }
  

  int OrderParamProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    // not implemented
    
    return -1;
  }

  //---PseudoRotation Class------------------------------------

  PseudoRotationProperty::PseudoRotationProperty(gcore::System &sys) :
    Property(sys)
  {
    d_type = "PseudoRotation";
    REQUIREDARGUMENTS = 1;
    d_sin36sin72 = sin(M_PI/5.0) + sin(2.0*M_PI/5.0);
    d_average = 0;
    d_zrmsd = 0;
    d_count = 0;
  }
  
  PseudoRotationProperty::~PseudoRotationProperty()
  {
  }
  
  void PseudoRotationProperty::parse(int count, std::string arguments[])
  {
    Property::parse(count, arguments);
    
    // it's a pseudo rotation, therefore 5 atoms needed
    if (d_atom.size() != 5)
      throw PseudoRotationProperty::Exception(
	      " wrong number of atoms for pseudo rotation.\n");
  }
  double PseudoRotationProperty::calc()
  {
    // first calculate the four dihedrals
    double torsion[5];
    torsion[0] = _calcDihedral(0,1,2,3);
    torsion[1] = _calcDihedral(1,2,3,4);
    torsion[2] = _calcDihedral(2,3,4,0);
    torsion[3] = _calcDihedral(3,4,0,1);
    torsion[4] = _calcDihedral(4,0,1,2);

    for(int i=0; i<5; ++i)
      if(torsion[i] > 180.0) torsion[i]-=360.0;
    
    double factor= (torsion[2] + torsion[4]) - (torsion[1] + torsion[3]);
    
    factor /= (2.0 * torsion[0] * d_sin36sin72);

    d_value = atan(factor)*180/M_PI;

    if(torsion[0] < 0.0) d_value += 180;

    d_average += d_value;
    d_zrmsd += pow(d_value-d_zvalue, 2);
    
    ++d_count;
    return d_value;

  }
  
  double PseudoRotationProperty::
  _calcDihedral(int const a, int const b, int const c, int const d)
  {
    gmath::Vec tmpA = atoms().pos(a) - atoms().pos(b);
    gmath::Vec tmpB = atoms().pos(d) - atoms().pos(c);
    gmath::Vec tmpC = atoms().pos(c) - atoms().pos(b);

    gmath::Vec p1 = tmpA.cross(tmpC);
    gmath::Vec p2 = tmpB.cross(tmpC);

    double cosphi = ((p1.dot(p2))/(p1.abs()*p2.abs()));

    double value = acos(cosphi)*180/M_PI;     

    gmath::Vec p3 = p1.cross(p2);
    if (p3.dot(tmpC)<0)
      value = 360 - value;   
    
    return value;
  }

  std::string PseudoRotationProperty::average()
  {
    std::ostringstream os;
    double z = sqrt(d_zrmsd / d_count);
    
    os << std::setw(20) << d_average / d_count
       << std::setw(20) << z;

    return os.str();
  }
  
  int PseudoRotationProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    // maybe give the residue number?
    if(mol_topo.resNum(atoms().atom(0)) == mol_topo.resNum(atoms().atom(1)) &&
       mol_topo.resNum(atoms().atom(0)) == mol_topo.resNum(atoms().atom(2)) &&
       mol_topo.resNum(atoms().atom(0)) == mol_topo.resNum(atoms().atom(3)) &&
       mol_topo.resNum(atoms().atom(0)) == mol_topo.resNum(atoms().atom(4)) )
      return mol_topo.resNum(atoms().atom(0));
    else
      return -1;
  }

  //---PuckerAmplitude Class------------------------------------

  PuckerAmplitudeProperty::PuckerAmplitudeProperty(gcore::System &sys) :
    PseudoRotationProperty(sys)
  {
    d_type = "PuckerAmplitude";
    //REQUIREDARGUMENTS = 1;
    //d_sin36sin72 = sin(M_PI/5.0) + sin(2.0*M_PI/5.0);
    //d_average = 0;
    //d_zrmsd = 0;
    //d_count = 0;
  }
  
  PuckerAmplitudeProperty::~PuckerAmplitudeProperty()
  {
  }
  
  double PuckerAmplitudeProperty::calc()
  {
    // first calculate the four dihedrals
    double torsion[5];
    torsion[0] = _calcDihedral(0,1,2,3);
    torsion[1] = _calcDihedral(1,2,3,4);
    torsion[2] = _calcDihedral(2,3,4,0);
    torsion[3] = _calcDihedral(3,4,0,1);
    torsion[4] = _calcDihedral(4,0,1,2);
    for(int i=0; i<5; ++i)
      if(torsion[i] > 180.0) torsion[i]-=360.0;
    
    double factor= (torsion[2] + torsion[4]) - (torsion[1] + torsion[3]);
    factor /= (2.0 * torsion[0] * d_sin36sin72);
    
    double pr = atan(factor);
    double t0 = torsion[0];
    
    if(torsion[0] < 0 ) {
      pr += M_PI;
    }
    
    d_value = t0 / cos(pr);
    
    d_average += d_value;
    d_zrmsd += pow(d_value-d_zvalue, 2);
    
    ++d_count;
    return d_value;

  }
  
} // utils


