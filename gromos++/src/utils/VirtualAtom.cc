#include <cassert>
#include <vector>
#include <sstream>
#include <iomanip>
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "Neighbours.h"
#include "VirtualAtom.h"


using utils::VirtualAtom;
using utils::VirtualAtom_i;
using gmath::Vec;
using namespace gcore;

static const double TETHCO=0.577350;
static const double TETHSI=0.816497;

class VirtualAtom_i{
  friend class VirtualAtom;
  
  // System that is restraint...
  const System &d_sys;
  // configuration as in table 2.6.4.1
  int d_config[4];
  // what molecule is restrained
  int d_mol;
  // type ICDR
  int d_type;
  int d_orient;
  
  static double DISH;
  static double DISC;
  

  VirtualAtom_i(const System &sys, int mol, int atom, int type, int orient): 
    d_sys(sys), d_mol(mol),d_type(type), d_orient(orient){
    d_config[0]=atom;
    
    for(int j=1;j<4;++j) d_config[j]=-1;
  }

  VirtualAtom_i(const System &sys, int mol, int atom, int type, int config[],  double dish, double disc, int orient): 
   d_sys(sys), d_mol(mol),d_type(type), d_orient(orient) {
   d_config[0]=atom;
    
   for(int j=1;j<4;++j) d_config[j]= config[j];
    DISH = dish;
    DISC = disc;

  }

  VirtualAtom_i(const VirtualAtom_i &v):
    d_sys(v.d_sys), d_mol(v.d_mol), d_type(v.d_type), d_orient(v.d_orient){

    for(int i=0; i<4; ++i) d_config[i] = v.d_config[i];
  }
  

  ~VirtualAtom_i(){}
};

double VirtualAtom_i::DISH = 0.1;
double VirtualAtom_i::DISC = 0.153;



VirtualAtom &VirtualAtom::operator=(const VirtualAtom &va){
   if (this != &va) {
     VirtualAtom_i tmp(*va.d_this);
     for(int i=0; i<4; ++i)  d_this->d_config[i] = tmp.d_config[i];
      d_this->d_mol = tmp.d_mol;
      d_this->d_type = tmp.d_type;
      d_this->d_orient = tmp.d_orient;
      d_this->DISH = tmp.DISH;
      d_this->DISC = tmp.DISC;
   }
 
   return *this;
}

VirtualAtom::VirtualAtom(const VirtualAtom &v) {
  if (this != &v)  d_this = new VirtualAtom_i(*v.d_this);
}


VirtualAtom::VirtualAtom(const System &sys, int mol, 
		 int atom, int type, int config[], double dish, double disc, int orientation)
  :  d_this(new VirtualAtom_i(sys, mol, atom, type, config, dish, disc, orientation)) 
{ }

VirtualAtom::VirtualAtom(const System &sys, int mol, 
			 int atom, int type, int orientation)
  :  d_this(new VirtualAtom_i(sys, mol, atom, type, orientation)) 
{
  Neighbours neigh(d_this->d_sys.mol(mol),atom);
  
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
    copy(neigh.begin(),neigh.end(),&(d_this->d_config[1]));
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
    copy(neigh.begin(),neigh.end(),&(d_this->d_config[1]));
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
    copy(neigh.begin(),neigh.end(),&(d_this->d_config[1]));
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
    copy(neigh.begin(),neigh.end(),&(d_this->d_config[1]));
    if(orientation==0 && d_this->d_config[1] > d_this->d_config[2])
      std::swap(d_this->d_config[1],d_this->d_config[2]);
    if(orientation==1 && d_this->d_config[1] < d_this->d_config[2])
      std::swap(d_this->d_config[1],d_this->d_config[2]);
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
    d_this->d_config[1]=neigh[0];   
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
    for(int l=0,m=1;l<3;l++)
      if(Neighbours(d_this->d_sys.mol(mol),neigh[l]).size()==1)
	d_this->d_config[m++]=neigh[l];

    break;
    
  default:
      
    std::ostringstream ss;
    ss << "Type " << type << " is unknown as a type of a virtual atom.";
    
    throw Exception(ss.str());

  }

}

VirtualAtom::~VirtualAtom(){
  delete d_this;
}


Vec VirtualAtom::pos()const
{
  Vec s,t;

  const Molecule &mol=d_this->d_sys.mol(d_this->d_mol);
  int *conf = d_this->d_config;

  const double &DISH=VirtualAtom_i::DISH;
  const double &DISC=VirtualAtom_i::DISC;
  
  
  switch(d_this->d_type){
    
  case 0: // explicit atom
  case 7: // rotating ring
    return mol.pos(d_this->d_config[0]);
    break;
       
  case 1: // CH1
    
    s=3.0 *mol.pos(conf[0])  - mol.pos(conf[1])
      - mol.pos(conf[2]) - mol.pos(conf[3]);
    return mol.pos(conf[0])+DISH/s.abs()*s;
    break;

  case 2: // aromatic H
       
    s=2.0*mol.pos(conf[0])-mol.pos(conf[1])
      -mol.pos(conf[2]);
    return mol.pos(conf[0])+DISH/s.abs()*s;
    break;
    
  case 3: // non-stereospecific CH2
    s=2.0*mol.pos(conf[0])-mol.pos(conf[1])
      -mol.pos(conf[2]);
    return mol.pos(conf[0])+DISH*TETHCO/s.abs()*s;
    break;
       
  case 4: // stereospecific CH2
    
    s = 2.0*mol.pos(conf[0])-mol.pos(conf[1]) -mol.pos(conf[2]);
    t = (mol.pos(conf[0])-mol.pos(conf[1])).cross(mol.pos(conf[0])-mol.pos(conf[2]));
    return mol.pos(conf[0]) + DISH*TETHCO/s.abs()*s + DISH * TETHSI / t.abs() * t;
    break;
       
  case 5: // CH3
       
    s = mol.pos(conf[0])-mol.pos(conf[1]);
    return mol.pos(conf[0]) + DISH / (3*s.abs()) * s;
    break;
       
  case 6: // non-stereospecific CH3 (Leu, Val)

    s = 2.0*mol.pos(conf[0])-mol.pos(conf[1])
      -mol.pos(conf[2]);
    return mol.pos(conf[0]) - TETHCO *(DISC+DISH/3.0)/s.abs()*s;
    break;

  case 8: // NH2-group (one pseudosite)

    s = 2.0*mol.pos(conf[0])-mol.pos(conf[1]) -mol.pos(conf[2]);
    return mol.pos(conf[0]) - (DISH * 0.5) * s/s.abs();
    break;

  default:
    throw Exception("Type code for virtual atom is not valid.");
    
  }

  return Vec(0,0,0);
  
}


void VirtualAtom::setDish(double dish)
{
  VirtualAtom_i::DISH=dish;
}
void VirtualAtom::setDisc(double disc)
{
  VirtualAtom_i::DISC=disc;
}

  
int VirtualAtom::type()const
{
  return d_this->d_type;
}

int VirtualAtom::mol()const{
  return d_this->d_mol;
}

int VirtualAtom::operator[](int i)const {
  assert(i<4&&i>=0);
  return d_this->d_config[i];
}


int VirtualAtom::orientation()const{
  return d_this->d_orient;
}
