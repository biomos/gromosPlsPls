// 	$Id$	
  
//---Property Class-----------------------------------

#include <cassert>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <typeinfo>

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
#include "../bound/Boundary.h"

#include "../gmath/Stat.h"
#include "Value.h"
#include "ExpressionParser.h"
#include "Property.h"

using namespace gcore;
using namespace std;
using namespace utils;

namespace utils
{
   Property::Property(gcore::System &sys, bound::Boundary * pbc)
     : d_type("Property"),
       d_atom(sys),
       d_sys(&sys),
       d_pbc(pbc)
   {
   }

  Property::~Property()
  {
  }

  void Property::parse(std::vector<std::string> const & arguments, int x)
  {
    if (int(arguments.size()) < REQUIREDARGUMENTS)
      throw Exception(" too few arguments.\n");

    if (REQUIREDARGUMENTS > 0)
      parseAtoms(arguments[0], x);

    d_arg.resize(arguments.size() - 1, Value(0.0));

    for(unsigned int i=1; i<arguments.size(); ++i)
      d_arg[i-1].parse(arguments[i]);
  }
  
  void Property::parse(AtomSpecifier const & atmspc)
  {
    d_atom=atmspc;
  }

  void Property::parseAtoms(std::string atoms, int x)
  {
    d_atom.addSpecifier(atoms, x);
  }

  std::string Property::toTitle()const
  { 
    std::ostringstream os;
    os << d_type << "%" << d_atom.toString()[0];
    return os.str();
  }

  std::ostream &operator<<(std::ostream &os, Property const & p)
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

  std::string Property::average()const
  {
    std::ostringstream os;
    os.precision(8);
    os.setf(std::ios::fixed, std::ios::floatfield);

    if (d_scalar_stat.n()){
      os << std::setw(15) << d_scalar_stat.ave()
	 << std::setw(15) << d_scalar_stat.rmsd()
	 << std::setw(15) << d_scalar_stat.ee();
    }
    if (d_vector_stat.n()){
      // std::cerr << "getting vector stat: ave" << std::endl;
      os << std::setw(15) << gmath::v2s(d_vector_stat.ave());
      // std::cerr << "getting vector stat: rmsd" << std::endl;      
      os << std::setw(15) << gmath::v2s(d_vector_stat.rmsd());
      // std::cerr << "getting vector stat: ee" << std::endl;      
      // os << std::setw(15) << gmath::v2s(d_vector_stat.ee());
    }
    
    return os.str();
  }

  void Property::addValue(Value const & v)
  {
    switch(v.type()){
      case val_scalar:
	d_scalar_stat.addval(v.scalar());
	break;
      case val_vector:
      case val_vecspec:
	d_vector_stat.addval(v.vec());
	break;
      default:
	throw Exception("bad value type");
    }
  }

  //---AverageProperty Class------------------------------------

  AverageProperty::AverageProperty
  (
   gcore::System &sys,
   bound::Boundary * pbc)
    : Property(sys, pbc)
  {
    d_type = "Average";
    REQUIREDARGUMENTS = 0;
  }

  Value const & AverageProperty::calc()
  {
    // empty
    d_single_scalar_stat = gmath::Stat<double>();
    d_single_vector_stat = gmath::Stat<gmath::Vec>();
    
    for(unsigned int i=0; i<d_property.size(); ++i){
      Value const & v = d_property[i]->calc();

      switch(v.type()){
	case val_scalar:
	  d_single_scalar_stat.addval(v.scalar());
	  break;
	case val_vector:
	case val_vecspec:
	  d_single_vector_stat.addval(v.vec());
	  break;
	default:
	  throw Exception("wrong value type");
      }
    }
    
    if (d_single_scalar_stat.n())
      d_scalar_stat.addval(d_single_scalar_stat.ave());
    if (d_single_vector_stat.n())
      d_vector_stat.addval(d_single_vector_stat.ave());

    if (d_single_scalar_stat.n() < d_single_vector_stat.n())
      d_value = d_single_scalar_stat.ave();
    else
      d_value = d_single_vector_stat.ave();
    
    return d_value;
  }

  std::string AverageProperty::toTitle()const
  {
    std::ostringstream os;
    os << "<"; 
    bool first = true;
    for(unsigned int i=0; i<d_property.size(); ++i){
      if (!first)
	os << ",";
      first = false;
      os << d_property[i]->toTitle();
    }
    os << ">";
    return os.str();
  }

  std::string AverageProperty::toString()const
  {
    std::ostringstream os;
    
    if (d_single_scalar_stat.n()){
      
      Value av = Value(d_single_scalar_stat.ave());
      Value rmsd = Value(d_single_scalar_stat.rmsd());
      Value ee = Value(d_single_scalar_stat.ee());
      
      os << av.toString() << " " << rmsd.toString() << " " << ee.toString();
      
      if (d_single_vector_stat.n()) os << "\t";
    }

    if (d_single_vector_stat.n()){
      
      Value av = Value(d_single_vector_stat.ave());
      Value rmsd = Value(d_single_vector_stat.rmsd());
      Value ee = Value(d_single_vector_stat.ee());
      
      os << av.toString() << " " << rmsd.toString() << " " << ee.toString();
    }
    return os.str();
  }
    
  //---DistanceProperty Class------------------------------------

  DistanceProperty::DistanceProperty
  (
   gcore::System &sys,
   bound::Boundary * pbc)
    : Property(sys, pbc)
  {
    d_type = "Distance";
    REQUIREDARGUMENTS = 1;
  }
  
  void DistanceProperty::parse(std::vector<std::string> const & arguments, int x)
  {
    Property::parse(arguments, x);

    // for a distance, we should just have to atoms here...
    if (atoms().size() != 2)
      throw Exception("wrong number of atoms for a distance.\n");
  }

  void DistanceProperty::parse(AtomSpecifier const & atmspc)
  {
    // for a distance, we should just have to atoms here...
    if (atmspc.size() != 2)
      throw Exception("wrong number of atoms for a distance.\n");
    Property::parse(atmspc);
  }

  Value const & DistanceProperty::calc()
  {
    gmath::Vec tmp = atoms().pos(0) - 
      d_pbc->nearestImage(atoms().pos(0), atoms().pos(1), d_sys->box());

    const double d = tmp.abs();
    
    d_value = d;
    addValue(d_value);

    return d_value;
  }

  int DistanceProperty::findTopologyType
  (
   gcore::MoleculeTopology const &mol_topo
   )
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

  AngleProperty::AngleProperty(gcore::System &sys, bound::Boundary * pbc) 
    : Property(sys, pbc)
  {
    d_type = "Angle";
    REQUIREDARGUMENTS = 1;
  }
  
  void AngleProperty::parse(std::vector<std::string> const & arguments, int x)
  {
    Property::parse(arguments, x);
    
    // it's an angle, therefore 3 atoms
    if (atoms().size() != 3)
      throw Exception("wrong number of atoms for an angle.\n");
  }

  void AngleProperty::parse(AtomSpecifier const & atmspc)
  {
    // it's an angle, therefore 3 atoms
    if (atmspc.size() != 3)
      throw Exception("wrong number of atoms for an angle.\n");
    Property::parse(atmspc);
  }
  
  Value const & AngleProperty::calc()
  {
    gmath::Vec tmpA = atoms().pos(0) 
      - d_pbc->nearestImage(atoms().pos(0), atoms().pos(1), d_sys->box());
    gmath::Vec tmpB = atoms().pos(2)
      - d_pbc->nearestImage(atoms().pos(2), atoms().pos(1), d_sys->box());

    const double d = acos((tmpA.dot(tmpB))/(tmpA.abs()*tmpB.abs()))*180/M_PI;

    d_value = d;
    addValue(d_value);
    
    return d_value;
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

  TorsionProperty::TorsionProperty(gcore::System &sys, bound::Boundary * pbc) :
    Property(sys, pbc)
  {
    d_type = "Torsion";
    REQUIREDARGUMENTS = 1;
  }
  
  TorsionProperty::~TorsionProperty()
  {
  }
  
  void TorsionProperty::parse(std::vector<std::string> const & arguments, int x)
  {
    Property::parse(arguments, x);
    
    // it's a torsion, therefore 4 atoms needed
    if (d_atom.size() != 4)
      throw Exception("wrong number of atoms for torsion.\n");
  }

  void TorsionProperty::parse(AtomSpecifier const & atmspc)
  {
    // it's a torsion, therefore 4 atoms needed
    if (atmspc.size() != 4)
      throw Exception("wrong number of atoms for torsion.\n");
    Property::parse(atmspc);
  }
  
  Value const & TorsionProperty::calc()
  {
    gmath::Vec tmpA = atoms().pos(0) - d_pbc->nearestImage(atoms().pos(0), atoms().pos(1), d_sys->box());
    gmath::Vec tmpB = atoms().pos(3) - d_pbc->nearestImage(atoms().pos(3), atoms().pos(2), d_sys->box());
    gmath::Vec tmpC = atoms().pos(2) - d_pbc->nearestImage(atoms().pos(2), atoms().pos(1), d_sys->box());

    gmath::Vec p1 = tmpA.cross(tmpC);
    gmath::Vec p2 = tmpB.cross(tmpC);

    double cosphi = ((p1.dot(p2))/(p1.abs()*p2.abs()));

    if(cosphi > 1.0) cosphi=1.0;
    if(cosphi <-1.0) cosphi=-1.0;
    double d = acos(cosphi)*180/M_PI;     

    gmath::Vec p3 = p1.cross(p2);
    if (p3.dot(tmpC)<0)
      d = 360 - d;   
    
    d_value = d;
    addValue(d_value);
    
    return d_value;
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

  OrderProperty::OrderProperty(gcore::System &sys, bound::Boundary * pbc) :
    Property(sys, pbc),
    d_axis(0)
  {
    // the first one we read ourselves...
    REQUIREDARGUMENTS = 1;
    d_type = "Order";
  }
  
  OrderProperty::~OrderProperty()
  {
  }
  
  void OrderProperty::parse(std::vector<std::string> const & arguments, int x)
  {
    if (arguments[0] == "x") d_axis = Vec(1,0,0);
    else if (arguments[0] == "y") d_axis = Vec(0,1,0);
    else if (arguments[0] == "z") d_axis = Vec(0,0,1);
    else
      throw Exception("axis specification wrong.\n");

    std::vector<std::string> my_args = arguments;
    my_args.erase(my_args.begin());
    
    Property::parse(my_args, x);
    
    if (d_atom.size() != 2)
      throw Exception("wrong number of atoms for an order property.\n");
  }
  
  Value const & OrderProperty::calc()
  {
    gmath::Vec tmpA = atoms().pos(0) - d_pbc->nearestImage(atoms().pos(0), atoms().pos(1), d_sys->box());

    double d = acos((tmpA.dot(d_axis))/(tmpA.abs()*d_axis.abs()))*180/M_PI;
    
    d_value = d;
    addValue(d_value);
    
    return d_value;
  }
  
  std::string OrderProperty::toTitle()const
  {
    std::ostringstream os;
    os << "o%";
    if (d_axis[0]) os << "x";
    else if (d_axis[1]) os << "y";
    else os << "z";
    
    os << "%" << atoms().toString()[0];
    
    return os.str();
  }
  
  int OrderProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    // not implemented
    return -1;
  }

  //---VectorOrderProperty Class------------------------------------------------------------

  VectorOrderProperty::VectorOrderProperty(gcore::System &sys, bound::Boundary * pbc) :
    Property(sys, pbc),
    d_vec1(sys, pbc),
    d_vec2(sys, pbc)
  {
    // the first one we read ourselves...
    REQUIREDARGUMENTS = 2;
    d_type = "VectorOrder";
  }
  
  VectorOrderProperty::~VectorOrderProperty()
  {
  }
  
  void VectorOrderProperty::parse(std::vector<std::string> const & arguments, int x)
  {
    d_vec1.setSpecifier(arguments[0]);
    d_vec2.setSpecifier(arguments[1]);

    // nothing left to parse ...
    // Property::parse(count - 2, &arguments[1]);
  }
  
  Value const & VectorOrderProperty::calc()
  {
    gmath::Vec tmpA = d_vec2();
    gmath::Vec d_axis = d_vec1();

    const double d = acos((tmpA.dot(d_axis))/(tmpA.abs()*d_axis.abs()))*180/M_PI;
    
    d_value = d;
    addValue(d_value);

    return d_value;
  }
  
  std::string VectorOrderProperty::toTitle()const
  {

    std::ostringstream s;
    s << "vo%"
      << d_vec1.toString() << "%" << d_vec2.toString();
    
    return s.str();
  }
  
  int VectorOrderProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    // not implemented
    return -1;
  }

  //---OrderParamProperty Class------------------------------------------------------------

  OrderParamProperty::OrderParamProperty(gcore::System &sys, bound::Boundary * pbc) :
    Property(sys, pbc),
    d_axis(0)
  {
    // the first one we read ourselves...
    REQUIREDARGUMENTS = 1;
    d_type = "OrderParam";
  }
  
  OrderParamProperty::~OrderParamProperty()
  {
  }
  
  void OrderParamProperty::parse(std::vector<std::string> const & arguments, int x)
  {
    if (arguments[0] == "x") d_axis = Vec(1,0,0);
    else if (arguments[0] == "y") d_axis = Vec(0,1,0);
    else if (arguments[0] == "z") d_axis = Vec(0,0,1);
    else
      throw Exception("axis specification wrong.\n");

    std::vector<std::string> my_args = arguments;
    my_args.erase(my_args.begin());
    Property::parse(my_args, x);
    
    if (d_atom.size() != 2)
      throw Exception("wrong number of atoms for an order property.\n");
  }
  
  Value const & OrderParamProperty::calc()
  {
    gmath::Vec tmpA = atoms().pos(0) - d_pbc->nearestImage(atoms().pos(0), atoms().pos(1), d_sys->box());

    const double cosa = tmpA.dot(d_axis)/(tmpA.abs()*d_axis.abs());
    const double d = 0.5 * (3 * cosa * cosa - 1);
    
    d_value = d;
    addValue(d_value);
    
    return d_value;
  }
  
  std::string OrderParamProperty::toTitle()const
  {

    std::ostringstream os;
    os << "op%";
    if (d_axis[0]) os << "x";
    else if (d_axis[1]) os << "y";
    else os << "z";
    
    os << "%" << atoms().toString()[0];
    return os.str();
  }

  int OrderParamProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    // not implemented
    return -1;
  }

  //---VectorOrderParamProperty Class------------------------------------------------------------

  VectorOrderParamProperty::VectorOrderParamProperty(gcore::System &sys, bound::Boundary * pbc) :
    Property(sys, pbc),
    d_vec1(sys, pbc),
    d_vec2(sys, pbc)
  {
    // the first one we read ourselves...
    REQUIREDARGUMENTS = 2;
    d_type = "VectorOrderParam";
  }
  
  VectorOrderParamProperty::~VectorOrderParamProperty()
  {
  }
  
  void VectorOrderParamProperty::parse(std::vector<std::string> const & arguments, int x)
  {
    d_vec1.setSpecifier(arguments[0]);
    d_vec2.setSpecifier(arguments[1]);

    // nothing left to parse ...
    // Property::parse(count - 1, &arguments[1]);
  }
  
  Value const & VectorOrderParamProperty::calc()
  {
    gmath::Vec tmpA = d_vec2();
    gmath::Vec d_axis = d_vec1();

    const double cosa = tmpA.dot(d_axis)/(tmpA.abs()*d_axis.abs());
    const double d = 0.5 * (3 * cosa * cosa - 1);
    
    d_value = d;
    addValue(d_value);

    return d_value;
  }
  
  std::string VectorOrderParamProperty::toTitle()const
  {

    std::ostringstream s;
    s << "vop%"
      << d_vec1.toString() << "%" << d_vec2.toString();
    
    return s.str();
  }

  int VectorOrderParamProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    // not implemented
    return -1;
  }

  //---PseudoRotation Class------------------------------------

  PseudoRotationProperty::PseudoRotationProperty(gcore::System &sys, bound::Boundary * pbc) :
    Property(sys, pbc),
    d_sin36sin72(::sin(M_PI/5.0) + ::sin(2.0*M_PI/5.0))
  {
    d_type = "PseudoRotation";
    REQUIREDARGUMENTS = 1;
  }
  
  PseudoRotationProperty::~PseudoRotationProperty()
  {
  }
  
  void PseudoRotationProperty::parse(std::vector<std::string> const & arguments, int x)
  {
    Property::parse(arguments, x);
    
    // it's a pseudo rotation, therefore 5 atoms needed
    if (d_atom.size() != 5)
      throw Exception(
	      " wrong number of atoms for pseudo rotation.\n");
  }

  Value const & PseudoRotationProperty::calc()
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

    double d = atan(factor)*180/M_PI;

    if(torsion[0] < 0.0) d += 180;

    d_value = d;
    addValue(d_value);
    
    return d_value;
  }
  
  double PseudoRotationProperty::
  _calcDihedral(int const a, int const b, int const c, int const d)
  {
    gmath::Vec tmpA = atoms().pos(a) - d_pbc->nearestImage(atoms().pos(a), atoms().pos(b), d_sys->box());
    gmath::Vec tmpB = atoms().pos(d) - d_pbc->nearestImage(atoms().pos(d), atoms().pos(c), d_sys->box());
    gmath::Vec tmpC = atoms().pos(c) - d_pbc->nearestImage(atoms().pos(c), atoms().pos(b), d_sys->box());

    gmath::Vec p1 = tmpA.cross(tmpC);
    gmath::Vec p2 = tmpB.cross(tmpC);

    double cosphi = ((p1.dot(p2))/(p1.abs()*p2.abs()));

    double value = acos(cosphi)*180/M_PI;     

    gmath::Vec p3 = p1.cross(p2);
    if (p3.dot(tmpC)<0)
      value = 360 - value;   
    
    return value;
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

  PuckerAmplitudeProperty::PuckerAmplitudeProperty(gcore::System &sys, bound::Boundary * pbc) :
    PseudoRotationProperty(sys, pbc)
  {
    d_type = "PuckerAmplitude";
    //REQUIREDARGUMENTS = 1;
    //d_sin36sin72 = sin(M_PI/5.0) + sin(2.0*M_PI/5.0);
  }
  
  PuckerAmplitudeProperty::~PuckerAmplitudeProperty()
  {
  }
  
  Value const & PuckerAmplitudeProperty::calc()
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
    
    const double d = t0 / ::cos(pr);
    
    d_value = d;
    addValue(d_value);
    
    return d_value;
  }

  //---ExpressionProperty Class------------------------------------------------------------

  ExpressionProperty::ExpressionProperty(gcore::System &sys, bound::Boundary * pbc) :
    Property(sys, pbc),
    d_parser(&sys, pbc)
  {
    // the first one we read ourselves...
    REQUIREDARGUMENTS = 1;
    d_type = "Expression";
  }
  
  ExpressionProperty::~ExpressionProperty()
  {
  }
  
  void ExpressionProperty::parse(std::vector<std::string> const & arguments, int x)
  {
    d_expr.clear();
    d_var["x"] = x;
    d_parser.parse_expression(arguments[0], d_var, d_expr);

    // nothing left to parse ...
    // Property::parse(count - 2, &arguments[1]);
  }
  
  Value const & ExpressionProperty::calc()
  {
    d_value = d_parser.calculate(d_expr, d_var);
    addValue(d_value);

    return d_value;
  }
  
  std::string ExpressionProperty::toTitle()const
  {
    std::ostringstream s;
    s << "expr";
    return s.str();
  }
  
  int ExpressionProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    // not implemented
    return -1;
  }

//---HBProperty Class------------------------------------------------------------------
  
  HBProperty::HBProperty(gcore::System &sys, bound::Boundary * pbc) :
    Property(sys, pbc),
    /**  HB distances (2 if 3 center)
     */
    d1_hb(sys, pbc), d2_hb(sys, pbc),
    /** HB angles (3 if 3 center)
     */
    a1_hb(sys, pbc), a2_hb(sys, pbc), a3_hb(sys, pbc),
    /** HB improper (1 if 3 center)
     */
    i1_hb(sys, pbc),
    /** auxiliary atom specifier
     */
    as(sys)
  {
    d_type = "HB";
    REQUIREDARGUMENTS = 1;
  }
  
  HBProperty::~HBProperty()
  {
  }
  
  void HBProperty::parse(std::vector<std::string> const & arguments, int x)
  {
    if (int(arguments.size()) < REQUIREDARGUMENTS)
      throw Exception(" too few arguments.\n");
    
    // 'as' order should be D-H-A1(-A2)
    as.clear();
    as.addSpecifier(arguments[0]);
    
    if (as.mass(1) != 1.00800)
      throw Exception(" the second atom must be a Hydrogen.\n");
    
    if (as.size() == 3) {   // no 3 center hb
      is3c=false;
      a1_hb.parse(as);
      as.removeAtom(0);
      d1_hb.parse(as);

      d_arg.resize(2, Value(0.0));
     
      if(arguments.size() == 1) {
        d_arg[0].parse("0.25");   //default values for hbond
        d_arg[1].parse("135");    // max dist && min ang
      }
      else
        {
          for(unsigned int i=1; i<arguments.size(); ++i)
            d_arg[i-1].parse(arguments[i]); 
        }
    }
    else if (as.size() == 4) { // 3 center hb
      is3c=true;
      AtomSpecifier tmp; // need to reorder 'as'
      tmp.addAtom(as.mol(1), as.atom(1));
      tmp.addAtom(as.mol(3), as.atom(3));
      tmp.addAtom(as.mol(2), as.atom(2));
      tmp.addAtom(as.mol(0), as.atom(0));
      i1_hb.parse(tmp);

      tmp.removeAtom(3); tmp.removeAtom(2);
      d2_hb.parse(tmp);

      tmp.clear();
      tmp.addAtom(as.mol(2), as.atom(2));
      tmp.addAtom(as.mol(1), as.atom(1));
      tmp.addAtom(as.mol(3), as.atom(3));
      a3_hb.parse(tmp);
      
      as.removeAtom(2);
      a2_hb.parse(as);
      
      as.clear();
      as.addSpecifier(arguments[0]);
      as.removeAtom(3);
      a1_hb.parse(as);
      
      as.removeAtom(0);
      d1_hb.parse(as);

      d_arg.resize(4, Value(0.0));

      if(arguments.size() == 1) {
        d_arg[0].parse("0.27");   //default values for hbond
        d_arg[1].parse("90");     // max dist && min ang
        d_arg[2].parse("340");    // & min summ & max imp
        d_arg[3].parse("15");
      }
      else
        {
          for(unsigned int i=1; i<arguments.size(); ++i)
            d_arg[i-1].parse(arguments[i]); 
        }  
    }
    else
      throw Exception(" wrong number of atoms for a HB.\n");
      
    as.clear();
    as.addSpecifier(arguments[0]);
    
  }

  Value const & HBProperty::calc()
  {
    if(is3c)
      if (d1_hb.calc().scalar() <= d_arg[0].scalar() &&
          d2_hb.calc().scalar() <= d_arg[0].scalar() &&
          a1_hb.calc().scalar() >= d_arg[1].scalar() &&
          a2_hb.calc().scalar() >= d_arg[1].scalar() &&
          (a1_hb.calc().scalar() + 
           a2_hb.calc().scalar() +
           a3_hb.calc().scalar()) >= d_arg[2].scalar() &&
          i1_hb.calc().scalar() <= d_arg[3].scalar())
        d_value=1;
      else
        d_value=0;  
    else   // no 3c HB
      if (d1_hb.calc().scalar() <= d_arg[0].scalar() && 
          a1_hb.calc().scalar() >= d_arg[1].scalar())
        d_value=1;
      else
        d_value=0;
    addValue(d_value);

    return d_value;
  }

  int HBProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo)
  {
    if(is3c) 
      return 1;
    else
      return 0;
  }  
  
  std::string HBProperty::toTitle()const
  {
    std::ostringstream os;
    if(is3c)
      os << d_type << "%" << as.toString()[0]
         << "%" << d_arg[0].scalar()
         << "%" << d_arg[1].scalar()
         << "%" << d_arg[2].scalar()
         << "%" << d_arg[3].scalar(); 
    else
      os << d_type << "%" << as.toString()[0]
         << "%" << d_arg[0].scalar()
         << "%" << d_arg[1].scalar(); 
    return os.str();
  }

} // utils


