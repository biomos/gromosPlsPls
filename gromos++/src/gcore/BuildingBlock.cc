// gcore_BuildingBlock.cc

#include <cassert>
#include <new>
#include <string>
#include <sstream>
#include "BuildingBlock.h"
#include "BbSolute.h"
#include "SolventTopology.h"
#include "../gromos/Exception.h"

using gcore::BuildingBlock;

BuildingBlock::BuildingBlock():
  d_bb(),
  d_be(),
  d_bs(),
  d_fpepsi(),
  d_hbar(),
  d_ffcode("_no_FORCEFIELD_block_given_"),
  d_linkExclusions(),
  d_empty(true)
{}


BuildingBlock::BuildingBlock(const BuildingBlock &bld):
  d_bb(bld.d_bb.size()),
  d_be(bld.d_be.size()),
  d_bs(bld.d_bs.size()),
  d_fpepsi(bld.Fpepsi()),
  d_hbar(bld.Hbar()),
  d_ffcode(bld.ForceField()),
  d_linkExclusions(bld.LinkExclusions()),
  d_empty(false)
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

void BuildingBlock::addBuildingBlock(const BuildingBlock &bld)
{
  // check force field code, fpepsi and hbar
  if(d_ffcode!="_no_FORCEFIELD_block_given_" && d_ffcode!=bld.ForceField())
    throw gromos::Exception("BuildingBlock", "Force-field code of building block files"
			    " are not identical\n" + d_ffcode + " and " +bld.ForceField());
  if(!d_empty){
    if(d_fpepsi!=bld.Fpepsi()){
      std::ostringstream os;
      os << "Value of FPEPSI is not identical in building block files\n"
	 << d_fpepsi << " != " << bld.Fpepsi() << std::endl;
      throw gromos::Exception("BuildingBlock", os.str());
    }
    if(d_hbar!=bld.Hbar()){
      std::ostringstream os;
      os << "Value of HBAR is not identical in building block files\n"
	 << d_hbar << " != " << bld.Hbar() << std::endl;
      throw gromos::Exception("BuildingBlock", os.str());
    }
  }
  else{
    d_fpepsi=bld.Fpepsi();
    d_hbar=bld.Hbar();
    d_ffcode=bld.ForceField();
    d_linkExclusions=bld.LinkExclusions();
    d_empty=false;
  }
  
  // now add all the individual building blocks
  for(int i=0; i< bld.numBbSolutes(); ++i)
    addBbSolute(bld.bb(i));
  for(int i=0; i< bld.numBbEnds(); ++i)
    addBbEnd(bld.be(i));
  for(int i=0; i< bld.numBbSolvents(); ++i)
    addBbSolvent(bld.bs(i));
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







