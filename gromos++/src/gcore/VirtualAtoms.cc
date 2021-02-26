// gcore_VirtualAtoms.cc
/**
 * Class VirtualAtoms
 */

#include <cassert>
#include <vector>
#include <set>
#include "../utils/VirtualAtom.h"
#include "../utils/AtomSpecifier.h"
#include "System.h"
#include "VirtualAtoms.h"

using namespace std;
using gcore::VirtualAtoms;

VirtualAtoms::VirtualAtoms(){}

VirtualAtoms::VirtualAtoms(utils::AtomSpecifier as)
{
  for(unsigned int i =0; i< as.size(); i++){
    if(as.atom()[i]->type()==utils::spec_virtual){ 
      utils::VirtualAtom va(*(as.sys()), 
 	                    as.atom()[i]->conf(), 
                            static_cast<utils::VirtualAtom::virtual_type>(as.atom()[i]->virtualType()));
      va.setDish(d_dish);
      va.setDisc(d_disc);
      d_vas.push_back(va);
      d_iac.push_back(-1);
      d_charge.push_back(0.0);
    } else if(as.atom()[i]->type()==utils::spec_solute){
      utils::VirtualAtom va(*(as.sys()),
                            as.mol(i),
                            as.atom(i),
                            static_cast<utils::VirtualAtom::virtual_type>(0));
      d_vas.push_back(va);
      d_iac.push_back(as.iac(i));
      d_charge.push_back(as.charge(i));
    } else {
      throw gromos::Exception("VirtualAtoms", "Don't know how to handle atomspecifier");
    }
  }
}
    
VirtualAtoms::VirtualAtoms(const VirtualAtoms &vas):
   d_vas(vas.d_vas),
   d_iac(vas.d_iac),
   d_charge(vas.d_charge),
   d_dish(vas.d_dish),
   d_disc(vas.d_disc)
{
}

VirtualAtoms::~VirtualAtoms()
{
  d_vas.erase(d_vas.begin(), d_vas.end());
}
void VirtualAtoms::addVirtualAtom(utils::AtomSpecifier as, int iac, double charge )
{
for(unsigned int i =0; i< as.size(); i++){
    if(as.atom()[i]->type()==utils::spec_virtual){
      utils::VirtualAtom va(*(as.sys()),
                            as.atom()[i]->conf(),
                            static_cast<utils::VirtualAtom::virtual_type>(as.atom()[i]->virtualType()));
      d_vas.push_back(va);
      d_iac.push_back(iac);
      d_charge.push_back(charge);
    } else if(as.atom()[i]->type()==utils::spec_solute){
      utils::VirtualAtom va(*(as.sys()),
                            as.mol(i),
                            as.atom(i),
                            static_cast<utils::VirtualAtom::virtual_type>(0));
      d_vas.push_back(va);
      d_iac.push_back(as.iac(i));
      d_charge.push_back(as.charge(i));
    } else {
      throw gromos::Exception("VirtualAtoms", "Don't know how to handle atomspecifier");
    }
  }
}


void VirtualAtoms::addVirtualAtom(System &sys, std::vector<int> conf, int type,
                                  double dish, double disc, 
                                  int iac, double charge )
{
  utils::VirtualAtom va(sys, static_cast<utils::VirtualAtom::virtual_type>(type), conf, dish, disc, 0);
  d_vas.push_back(va);
  d_iac.push_back(iac);
  d_charge.push_back(charge);
}
unsigned int VirtualAtoms::numVirtualAtoms()const{
     return d_vas.size();
}
utils::VirtualAtom VirtualAtoms::atom(int i){
    assert(i < d_vas.size());
    return d_vas[i];
}
const utils::VirtualAtom VirtualAtoms::atom(int i)const{
    assert(i < d_vas.size());
    return d_vas[i];
}
void VirtualAtoms::setDis(double dish, double disc){
  d_disc = disc;
  d_dish = dish;
}
void VirtualAtoms::setSystem(gcore::System &sys){
  for(unsigned int i=0; i< d_vas.size(); i++){
    d_vas[i].setSystem(sys);
  }
}
