// gcore_BbLink.cc

#include <cassert>
#include <string>
#include <set>
#include <vector>
#include <map>
#include <new>
#include "AtomTopology.h"
#include "Exclusion.h"
#include "Bond.h"
#include "Angle.h"
#include "Improper.h"
#include "Dihedral.h"
#include "LJException.h"
#include "MoleculeTopology.h"
#include "Exclusion.h"
#include "BbSolute.h"
#include "BbLink.h"

using namespace std;
using namespace gcore;

BbLink::BbLink(const BbLink& mt)
{
  for(int i=0; i<mt.numAtoms(); i++)
    MoleculeTopology::addAtom(mt.atom(i));
  
  for(int i=0; i<mt.numPexcl(); i++)
    addPexcl(mt.pexcl(i));
  
  BondIterator bi(mt);
  for(;bi;++bi)
    MoleculeTopology::addBond(bi());
  
  BondDipoleIterator bdi(mt);
  for(;bdi;++bdi)
    MoleculeTopology::addDipoleBond(bdi());
  
  AngleIterator ai(mt);
  for(;ai;++ai)
    MoleculeTopology::addAngle(ai());
  
  DihedralIterator di(mt);
  for(;di;++di)
    MoleculeTopology::addDihedral(di());
  
  ImproperIterator ii(mt);
  for(;ii;++ii)
    MoleculeTopology::addImproper(ii());
  
  LJExceptionIterator lji(mt);
  for(;lji;++lji)
    MoleculeTopology::addLJException(lji());
  
  setResName(mt.resName());
  setRep(mt.rep());

  for(int i=0; i < mt.numAtoms(); i++)
    setLinkRes(i, mt.linkRes(i));
}
// Methods

BbLink &BbLink::operator=(const BbLink &mt){
  if (this != &mt){
    this->BbLink::~BbLink();
    new(this) BbLink(mt);
  }
  return *this;
}

void BbLink::setLinkRes(const int a, const int i){
  if(d_linkres.size() <= a)
   d_linkres.resize(a+1); 
  d_linkres[a]=i;
}
  
int BbLink::linkRes(const int a)const{
  return d_linkres[a];
}
