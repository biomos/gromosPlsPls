#include <cassert>
#include <vector>
#include <sstream>
#include <iomanip>
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "Neighbours.h"
#include "VirtualAtom.h"
#include "AtomSpecifier.h"

using utils::VirtualAtom;
using utils::VirtualAtom_i;
using gmath::Vec;
using namespace gcore;

static const double TETHCO=0.577350;
static const double TETHSI=0.816497;


//==============================================================================
// IMPLEMENTATION CLASS
//==============================================================================

class VirtualAtom_i{
  friend class VirtualAtom;
  
  const System *d_sys;

  // configuration as in table 2.6.4.1
  AtomSpecifier d_config;

  // type ICDR
  VirtualAtom::virtual_type d_type;
  
  double d_dish;
  double d_disc;

  int d_orient;
  
  /**
   * Constructor
   * probably the atom should go out... (as it is in config[0] ?!)
   */
  VirtualAtom_i(System &sys,
		int mol, int atom, VirtualAtom::virtual_type type, 
		std::vector<int> const &config,
		double dish, double disc, int orient)
    : d_sys(&sys), d_config(sys), d_type(type),
      d_dish(dish), d_disc(disc), d_orient(orient)
  {

    d_config.addAtom(mol, atom);
    
    // copy the config into the AtomSpecifier
    for(unsigned int i=1; i<config.size(); ++i)
      d_config.addAtom(mol, config[i]);

  }

  /**
   * Constructor
   * simpler, new style constructor.
   */
  VirtualAtom_i(const System &sys,
		AtomSpecifier const & config,
		VirtualAtom::virtual_type type,
		double dish, double disc,
		int orientation)
    : d_sys(&sys), d_config(config), d_type(type),
      d_dish(dish), d_disc(disc), d_orient(orientation)
  {
  }

  /**
   * Copy constructor
   */
  VirtualAtom_i(const VirtualAtom_i &v)
    : d_sys(v.d_sys), d_config(v.d_config),
      d_type(v.d_type),
      d_dish(v.d_dish), d_disc(v.d_disc),
      d_orient(v.d_orient)
  {
  }
  
  /**
   * Destructor
   */
  ~VirtualAtom_i(){}

  /**
   * change the system
   */
  void setSystem(gcore::System &sys)
  {
    d_sys = &sys;
    d_config.setSystem(sys);
  }
  
  
  ///////////////////////////////////////////////////////////////////////////
  // OLD FUNCTIONALITY
  ///////////////////////////////////////////////////////////////////////////
  /**
   * Constructor
   * @deprecated
   */
  VirtualAtom_i(System &sys, int mol, int atom,
		VirtualAtom::virtual_type type, int orient): 
    d_sys(&sys), d_config(sys), d_type(type), d_orient(orient){

    d_dish = 0.1;
    d_disc = 0.153;
    
    d_config.addAtom(mol, atom);
    // copy the config into the AtomSpecifier
    /*
      for(int i=1; i<4; ++i)
      d_config.addAtom(mol, -1);
    */

  }

};

//==============================================================================
// END OF IMPLEMENTATION CLASS
//==============================================================================


/**
 * Constructor
 * this one is used by the new programs
 */
VirtualAtom::VirtualAtom(System &sys,
			 AtomSpecifier const &spec,
			 virtual_type type,
			 double dish, double disc, int orientation)
  :  d_this(new VirtualAtom_i(sys, spec, type, dish, disc, orientation)) 
{ }


/**
 * Constructor
 * this one is used by the new programs
 */
VirtualAtom::VirtualAtom(System &sys, int mol, 
			 int atom, virtual_type type,
			 std::vector<int> const &config,
			 double dish, double disc, int orientation)
  :  d_this(new VirtualAtom_i(sys, mol, atom, type, config, dish, disc, orientation)) 
{ }

/**
 * copy constructor
 */
VirtualAtom::VirtualAtom(const VirtualAtom &v) {
  if (this != &v)  d_this = new VirtualAtom_i(*v.d_this);
}

/**
 * operator=()
 */
VirtualAtom &VirtualAtom::operator=(const VirtualAtom &va){
  if (this != &va) {
    // this is not exception safe... (should use swap)
    // like this, there is no advantage in using tmp!!!
    VirtualAtom_i tmp(*va.d_this);

    d_this->d_sys = tmp.d_sys;
    d_this->d_config = tmp.d_config;
    d_this->d_type = tmp.d_type;
    d_this->d_orient = tmp.d_orient;
    d_this->d_dish = tmp.d_dish;
    d_this->d_disc = tmp.d_disc;
  }
 
  return *this;
}


VirtualAtom::~VirtualAtom(){
  delete d_this;
}

////////////////////////////////////////////////////////////////////////////////
// methods

Vec VirtualAtom::pos()const
{
  Vec s,t;

  AtomSpecifier & spec = d_this->d_config;

  // const Molecule &mol=d_this->d_sys.mol(d_this->d_mol);

  const double &DISH = d_this->d_dish;
  const double &DISC = d_this->d_disc;
  
  switch(d_this->d_type){
    
    case 0: // explicit atom
    case 7: // rotating ring
      return *spec.coord(0);
       
    case 1: // CH1
    
      s = 3.0 * spec.pos(0)  - spec.pos(1)
	- spec.pos(2) - spec.pos(3);
      return spec.pos(0) + DISH / s.abs() * s;

    case 2: // aromatic H
       
      s = 2.0 * spec.pos(0) - spec.pos(1) - spec.pos(2);
      return spec.pos(0) + DISH / s.abs() * s;
    
    case 3: // non-stereospecific CH2
      s= 2.0 * spec.pos(0) - spec.pos(1) - spec.pos(2);
      return spec.pos(0) + DISH * TETHCO / s.abs() * s;
       
    case 4: // stereospecific CH2
    
      s = 2.0 * spec.pos(0) - spec.pos(1) - spec.pos(2);
      t = (spec.pos(0) - spec.pos(1)).cross(spec.pos(0) - spec.pos(2));
      return spec.pos(0) + DISH * TETHCO / s.abs() * s + DISH * TETHSI / t.abs() * t;
       
    case 5: // CH3
       
      s = spec.pos(0) - spec.pos(1);
      return spec.pos(0) + DISH / (3 * s.abs()) * s;
       
    case 6: // non-stereospecific CH3 (Leu, Val)

      s = 2.0 * spec.pos(0) - spec.pos(1) - spec.pos(2);
      return spec.pos(0) - TETHCO * (DISC + DISH / 3.0) / s.abs() * s;

    case 8: // NH2-group (one pseudosite)

      s = 2.0 * spec.pos(0) - spec.pos(1) - spec.pos(2);
      return spec.pos(0) - (DISH * 0.5) * s / s.abs();

    case 9: // (CH3)3-group (one psuedosite)

      s = spec.pos(0) - spec.pos(1);
      t = ( DISC + DISH/3.0 );
      return spec.pos(0) + t / (3 * s.abs()) * s;
      
    case COM: // type 100
      {
	double m = 0;
	gmath::Vec v = 0;
	
	for(int i=0; i<spec.size(); ++i){
	  v += spec.mass(i) * spec.pos(i);
	  m += spec.mass(i);
	}
	return v / m;
      }

    case COG: // type 101
      {
	gmath::Vec v = 0;
	
	for(int i=0; i<spec.size(); ++i){
	  v += spec.pos(i);
	}
	return v / spec.size();
      }
      
    default:
      throw Exception("Type code for virtual atom is not valid.");
    
  }

  return Vec(0,0,0);
  
}

std::string VirtualAtom::toString()const
{
  std::ostringstream os;
  os << "va(" << d_this->d_type << ",";
  os << d_this->d_config.toString()[0] << ")";
  return os.str();
}

void VirtualAtom::setDish(double dish)
{
  d_this->d_dish = dish;
}
void VirtualAtom::setDisc(double disc)
{
  d_this->d_disc = disc;
}

  
VirtualAtom::virtual_type VirtualAtom::type()const
{
  return d_this->d_type;
}

utils::AtomSpecifier & VirtualAtom::conf(){
  return d_this->d_config;
}


int VirtualAtom::orientation()const{
  return d_this->d_orient;
}

void VirtualAtom::setSystem(gcore::System &sys)
{
  d_this->setSystem(sys);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
// OLD FUNCTIONALITY
// should be removed as the old programs are removed
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * old constructor
 * @deprecated
 */
VirtualAtom::VirtualAtom(System &sys, int mol, 
			 int atom, virtual_type type, int orientation)
  :  d_this(new VirtualAtom_i(sys, mol, atom, type, orientation)) 
{
  Neighbours neigh(d_this->d_sys->mol(mol),atom);
  //cout << "va: do we get here " << type << " " << orientation << endl;
  
  switch(type){
    case 0:
    case 7:
      break;
    case 1:
      // stereospecific CH
      if(neigh.size()!=3){
	//ostrstream ss;
	std::ostringstream ss;
	ss << "Specifying type 1 for atom " << atom
	   << " of molecule " << mol 
	   << " does not make sense!";
      
	throw Exception(ss.str());
      }

      for(unsigned int i=0; i < neigh.size(); ++i){
	d_this->d_config.addAtom(mol, neigh[i]);
      }
    
      break;
    case 2:
      // aromatic CH1
      if(neigh.size()!=2){
	// ostrstream ss;
	std::ostringstream ss;
	ss << "Specifying type " << type << " for atom " << atom
	   << " of molecule " << mol 
	   << " does not make sense!";
	throw Exception(ss.str());
      }

      for(unsigned int i=0; i < neigh.size(); ++i){
	d_this->d_config.addAtom(mol, neigh[i]);
      }
      break;

    case 3:
      // non-stereospecific CH2
      if(neigh.size()!=2){
	std::ostringstream ss;
	ss << "Specifying type " << type << " for atom " << atom
	   << " of molecule " << mol 
	   << " does not make sense!";
      
	throw Exception(ss.str());
      }

      for(unsigned int i=0; i < neigh.size(); ++i){
	d_this->d_config.addAtom(mol, neigh[i]);
      }
      break;
    
    case 4:
      // stereospecific CH2: Look also for the orientation flag...
      // I define for orientation == 0  the atoms i, j, k with j < k
      //          for orientation == 1  the atoms i, j, k with j > k
      if(neigh.size()!=2){
	std::ostringstream ss;
	ss << "Specifying type " << type << " for atom " << atom
	   << " of molecule " << mol 
	   << " does not make sense!";
      
	throw Exception(ss.str());
      }

      if(orientation==0 && neigh[0] > neigh[1])
	std::swap(neigh[0], neigh[1]);

      if(orientation==1 && neigh[0] < neigh[1])
	std::swap(neigh[0], neigh[1]);

      for(unsigned int i=0; i < neigh.size(); ++i){
	d_this->d_config.addAtom(mol, neigh[i]);
      }

      break;

    case 5:
      // non-stereospecific CH3
      if(neigh.size()!=1){
	std::ostringstream ss;
	ss << "Specifying type " << type << " for atom " << atom
	   << " of molecule " << mol 
	   << " does not make sense!";
      
	throw Exception(ss.str());
      }
      d_this->d_config.addAtom(mol, neigh[0]);   
      break;
      
    case 6:
      // non-stereospecific CH3 groups (Val, Leu)
      if(neigh.size()!=3){
	std::ostringstream ss;
	ss << "Specifying type " << type << " for atom " << atom
	   << " of molecule " << mol 
	   << " does not make sense!";
      
	throw Exception(ss.str());
      }
      for(int l=0;l<3;l++)
	if(Neighbours(d_this->d_sys->mol(mol),neigh[l]).size()==1)
	  d_this->d_config.addAtom(mol, neigh[l]);
    
      break;
    
    case 8:
      // non-stereospecific NH2 groups (Asn, Gln)
      if(neigh.size()!=3){
	std::ostringstream ss;
	ss << "Specifying type " << type << " for atom " << atom
	   << " of molecule " << mol 
	   << " does not make sense!";
	
	throw Exception(ss.str());
      }
      for(int l=0;l<3;l++)
	if(Neighbours(d_this->d_sys->mol(mol),neigh[l]).size()==1)
	  d_this->d_config.addAtom(mol, neigh[l]);
      
      break;

    case 9:
      // non-stereospecific (CH3)3
      if(neigh.size()!=4){
	std::ostringstream ss;
	ss << "Specifying type " << type << " for atom " << atom
	   << " of molecule " << mol 
	   << " does not make sense!";
      
	throw Exception(ss.str());
      }
      // now, how do we find the neighbour that is NOT one of the methyls?
      // it is the one that is does have more than one neighbour
      for(int l=0;l<4;l++)
	if(Neighbours(d_this->d_sys->mol(mol),neigh[l]).size()!=1)
	  d_this->d_config.addAtom(mol, neigh[l]);

      break;

    default:
      
      std::ostringstream ss;
      ss << "Type " << type << " is unknown as a type of a virtual atom.";
      
      throw Exception(ss.str());
      
  }

}
