// fit_PositionUtils.cc

#include "PositionUtils.h"
#include "Reference.h"
#include "../gmath/Vec.h"
#include "../gmath/Matrix.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"

#include <vector>
#include <string>

using namespace gmath;
using namespace gcore;
using fit::PositionUtils;

Vec PositionUtils::com(const gcore::System &sys){
  // calculate the center of mass of the molecule
  double totalMass=0;
  Vec cm;

  for(int m=0;m<sys.numMolecules();++m)
    for (int i=0;i < sys.mol(m).numAtoms(); i++) {
      cm += sys.mol(m).pos(i) * sys.mol(m).topology().atom(i).mass();
      totalMass += sys.mol(m).topology().atom(i).mass();
    }

  cm = (1.0/totalMass)*cm;
  
  return cm;
}

Vec PositionUtils::cog(const gcore::System &sys){
  
  // calculate the center of mass of the molecule
  Vec cm;
  int atoms=0;

  for(int m=0;m<sys.numMolecules();++m)
    for (int i=0;i < sys.mol(m).numAtoms(); i++) {
      cm = cm + sys.mol(m).pos(i);
      ++atoms;
    }

  cm = (1.0/double(atoms))*cm;
  
  return cm;
}

Vec PositionUtils::com(const System &sys, const Reference &ref){
  double totalMass=0;
  Vec cm;

  for(int m=0;m<sys.numMolecules();++m)
    for (int i=0;i < sys.mol(m).numAtoms(); i++) {
      cm += 
	ref.weight(m,i)
	* sys.mol(m).topology().atom(i).mass()
	* sys.mol(m).pos(i) ;

      totalMass += ref.weight(m,i) *
	sys.mol(m).topology().atom(i).mass();

    }

  cm = (1.0/totalMass)*cm;
  return cm;
}

Vec PositionUtils::cog(const System &sys, const Reference &ref){
  Vec cg;

  for(int m=0;m<sys.numMolecules();++m)
    for (int i=0;i < sys.mol(m).numAtoms(); i++) {

     cg += 
	ref.weight(m,i)
	* sys.mol(m).pos(i) ;
    }

  return cg;
}

void PositionUtils::translate(gcore::System *sys, const gmath::Vec &v){

  for(int m=0;m<sys->numMolecules();++m)
    for(int i=0;i<sys->mol(m).numAtoms();++i)
      sys->mol(m).pos(i)=sys->mol(m).pos(i)+v;

}

void PositionUtils::rotate(gcore::System *sys, const gmath::Matrix &mat){
  for(int j=0;j<sys->numMolecules();++j)
    for(int i=0;i<sys->mol(j).numAtoms();++i)
      sys->mol(j).pos(i)=mat*sys->mol(j).pos(i);
}

void PositionUtils::shiftToCom(gcore::System *sys){
  translate(sys,-com(*sys));
}

void PositionUtils::shiftToCom(gcore::System *sys, const Reference &ref){
  translate(sys, -com(*sys, ref));
}

void PositionUtils::shiftToCog(gcore::System *sys){
  translate(sys,-cog(*sys));
}

void PositionUtils::shiftToCog(gcore::System *sys, const Reference &ref){
  translate(sys, -cog(*sys, ref));
}

