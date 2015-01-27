// fit_TranslationalFit.cc

#include <cassert>
#include <set>
#include "TranslationalFit.h"
#include "Reference.h"
#include "PositionUtils.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"

using fit::TranslationalFit;
using fit::Reference;
using fit::PositionUtils;
using gcore::System;
using gmath::Vec;

// static Vec com(const System &sys, const Reference &ref);
// static Vec cog(const System &sys, const Reference &ref);


class fit::TranslationalFit_i{
  friend class fit::TranslationalFit;
  Reference *d_ref;
  TranslationalFit_i(Reference *ref)
  {
    d_ref=ref;
  }
  ~TranslationalFit_i(){}
};

TranslationalFit::TranslationalFit(Reference *ref, centre_enum centre):
  d_this(new TranslationalFit_i(ref))
{
  if (centre == fit::cog)
    PositionUtils::translate(&ref->sys(),PositionUtils::cog(ref->sys(),*ref));
  else
    PositionUtils::translate(&ref->sys(),PositionUtils::com(ref->sys(),*ref));
}

TranslationalFit::~TranslationalFit(){
  delete d_this;
}

void TranslationalFit::fit(gcore::System *sys)const{
  PositionUtils::translate(sys,PositionUtils::cog(*sys,*d_this->d_ref));
}

/*
static Vec com(const System &sys, const Reference &ref){
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


static Vec cog(const System &sys, const Reference &ref){
  Vec cg;

  for(int m=0;m<sys.numMolecules();++m)
    for (int i=0;i < sys.mol(m).numAtoms(); i++) {

     cg += 
	ref.weight(m,i)
	* sys.mol(m).pos(i) ;
    }

  return cg;
}
*/
