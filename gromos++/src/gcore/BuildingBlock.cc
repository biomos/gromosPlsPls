// gcore_BuildingBlock.cc

#include <cassert>
#include "BuildingBlock.h"
#include "BbSolute.h"
#include "SolventTopology.h"
#include <new>
#include <string>

using gcore::BuildingBlock;

BuildingBlock::BuildingBlock():
  d_bb(),
  d_be(),
  d_bs(),
  d_fpepsi(),
  d_hbar(),
  d_linkExclusions()
{}


BuildingBlock::BuildingBlock(const BuildingBlock &bld):
  d_bb(bld.d_bb.size()),
  d_be(bld.d_be.size()),
  d_bs(bld.d_bs.size()),
  d_fpepsi(),
  d_hbar(),
  d_linkExclusions()
{
  for (unsigned int i=0; i<d_bb.size();++i){
    d_bb[i]= new BbSolute(bld.bb(i));
  }
  for (unsigned int i=0; i<d_be.size();++i){
      d_be[i] = new BbSolute(bld.be(i));
  }
  for (unsigned int i=0; i<d_bs.size();++i){
      d_bs[i] = new SolventTopology(bld.bs(i));
  }
  d_fpepsi=bld.Fpepsi();
  d_hbar=bld.Hbar();
  d_linkExclusions=bld.LinkExclusions();
}

BuildingBlock::~BuildingBlock(){
  for (unsigned int i=0; i<d_bb.size();++i){
    delete d_bb[i];
  }
  for (unsigned int i=0; i<d_bs.size();++i){
    delete d_bs[i];
  }
  for (unsigned int i=0; i<d_be.size();++i){
    delete d_be[i];
  }
}

BuildingBlock &BuildingBlock::operator=(const BuildingBlock &bld){
  if(this != &bld){
    delete this;
    new(this) BuildingBlock(bld);
  }
  return *this;
}

void BuildingBlock::addBbSolute(const BbSolute &bb){
  d_bb.push_back(new BbSolute(bb));
}

void BuildingBlock::addBbSolvent(const SolventTopology &bs){
    d_bs.push_back(new SolventTopology(bs));
}

void BuildingBlock::addBbEnd(const BbSolute &be){
    d_be.push_back(new BbSolute(be));
}

int BuildingBlock::findBb(std::string s)
{
  //loop over all solute building blocks
  for(unsigned int i=0; i<d_bb.size(); ++i){
    if(d_bb[i]->resName()==s) return i+1;
  }
  //or, if not found over the end-building blocks
  for(unsigned int i=0; i<d_be.size(); ++i){
    if(d_be[i]->resName()==s) return -i-1;
  }
  return 0;
  
}
int BuildingBlock::findBs(std::string s)
{
  for(unsigned int i=0; i<d_bs.size(); ++i){
    if(d_bs[i]->solvName()==s) return i+1;
  }
  return 0;
}







