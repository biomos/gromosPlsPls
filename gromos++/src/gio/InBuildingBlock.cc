// gio_InBuildingBlock.cc

#include "InBuildingBlock.h"
#include "Ginstream.h"
#include "../gcore/BbSolute.h"
#include "../gcore/BbEnd.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Constraint.h"
#include "../gcore/Dihedral.h"
#include "../gcore/Improper.h"
#include "../gcore/Exclusion.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/BuildingBlock.h"

#include <map>
#include <deque>
#include <set>
#include <cmath>

using namespace gcore;
using gio::InBuildingBlock_i;
using gio::InBuildingBlock;

// Implementation class
class InBuildingBlock_i{
  friend class InBuildingBlock;
  Ginstream d_gin;
  BuildingBlock d_bld;
  string d_name;
  void init();
  InBuildingBlock_i (const char *name):
    d_gin(name), d_name(name){
    init();
  }
  ~InBuildingBlock_i(){
    d_gin.close();
  }
};

// Constructors

InBuildingBlock::InBuildingBlock(string name):
  d_this(new InBuildingBlock_i(name.c_str())){
}

InBuildingBlock::~InBuildingBlock(){
  delete d_this;
}

const string &InBuildingBlock::title()const{
  return d_this->d_gin.title();
}

const BuildingBlock &InBuildingBlock::building()const{
  return d_this->d_bld;
}

void InBuildingBlock_i::init(){

  if(!d_gin)
    throw InBuildingBlock::Exception("Could not open Building block file "+d_gin.name()+".");

  // generic variables
  double d[4];
  int i[5], num;
  string s;

  // Topphyscon block:
  if(! d_gin.check("TOPPHYSCON"))
    throw InBuildingBlock::Exception("Building block file "+d_gin.name()+" is corrupted:\nNo TOPPHYSCON block!");
  
  d_gin >> d[0] >> d[1];
  d_bld.setFpepsi(d[0]);
  d_bld.setHbar(d[1]);

  if(! d_gin.check())
    throw InBuildingBlock::Exception("Building block file "+d_gin.name()+" is corrupted:\nTOPPHYSCON block is not OK!");

  // Linkexclusions Block 
  if(! d_gin.check("LINKEXCLUSIONS"))
    throw InBuildingBlock::Exception("Building block file "+d_gin.name()+" is corrupted:\nNo LINKEXCLUSIONS block!");
  d_gin >> i[0];
  d_bld.setLinkExclusions(i[0]);

  if(! d_gin.check())
    throw InBuildingBlock::Exception("Building block file "+d_gin.name()+" is corrupted:\nLINKEXCLUSIONS block is not OK!");

  // MTBUILDBLSOLUTE blocks
  int cntn=0;
  if(!d_gin.check("MTBUILDBLSOLUTE"))
      throw InBuildingBlock::Exception("Building block file "+d_gin.name()+" is corrupted:\nNo MTBUILDBLSOLUTE");
  else cntn=1;

  while(cntn){
 
    BbSolute bb;
    int na, npe;
    d_gin >> s;

    bb.setResName(s);
    d_gin >> na >> npe;
    for(int j = 0; j<npe; j++){
	Exclusion e;
	d_gin >> i[0] >> i[1];
	for(int k=0;k<i[1];k++){
	    d_gin >> i[2];
	    e.insert(i[2]-1);
	}
	bb.addPexcl(e);
    }

    for(int j = 0; j<na; j++){
	AtomTopology at;
	d_gin >> i[0];
        d_gin >> s;
        at.setName(s);
        d_gin >> i[0];
        at.setIac(--i[0]);
        d_gin >> d[0] >> d[1];
	// WARNING: in the building block we use the mass code, 
	//          in the AtomTopology we usually have real masses
	at.setMass(d[0]-1);
        at.setCharge(d[1]);
        d_gin >> i[0];
        at.setChargeGroup(i[0]);
        // The trailing atoms do not have exclusions specified.
	// these are specified as the preceding exclusions of the 
	// next residue
        if(j<na - npe){
	    Exclusion e;
	    d_gin >> i[0];
	    for(int k=0;k<i[0];k++){
		d_gin >> i[1];
		e.insert(--i[1]);
	    }
	    at.setExclusion(e);
	}
        bb.addAtom(at);
    }
    // Bonds etc.

    // Bonds
    d_gin >> num;
    for (int j=0; j<num; j++){
      d_gin >> i[0] >> i[1] >> i[2];
      Bond bond(--i[0],--i[1]);
      bond.setType(--i[2]);
      bb.addBond(bond);
    }

    // Angles
    d_gin >> num;
    for (int j=0; j<num; j++){
	d_gin >> i[0] >> i[1] >> i[2] >> i[3];
	Angle angle(--i[0],--i[1],--i[2]);
	angle.setType(--i[3]);
	bb.addAngle(angle);
    }

    // Impropers
    d_gin >> num;
    for (int j=0; j<num; j++){
	d_gin >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
	Improper improper(--i[0],--i[1],--i[2],--i[3]);
	improper.setType(--i[4]);
	bb.addImproper(improper);
    }

    // Dihedrals
    d_gin >> num;
    for (int j=0; j<num; j++){
	d_gin >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
	Dihedral dihedral(--i[0],--i[1],--i[2],--i[3]);
	dihedral.setType(--i[4]);
	bb.addDihedral(dihedral);
    }


    if(! d_gin.check())
      throw InBuildingBlock::Exception("Building block file "+d_gin.name()+" is corrupted:\n MTBUILDBLSOLUTE block "+bb.resName()+" is not OK!");
    d_bld.addBbSolute(bb);
    d_gin >> s;
    if(s=="MTBUILDBLSOLUTE") cntn=1;
    else cntn=0; 
  }

  if(s=="MTBUILDBLSOLVENT") cntn=1;
  // MTBUILDBLSOLVENT blocks
  while(cntn){
 
    SolventTopology bs;
    int na;
    d_gin >> s;
    bs.setSolvName(s);
    d_gin >> na;

    for(int j = 0; j<na; j++){
	AtomTopology at;
	d_gin >> i[0];
        d_gin >> s;
        at.setName(s);
        d_gin >> i[0];
        at.setIac(--i[0]);
        d_gin >> d[0] >> d[1];
	// WARNING: in the building block we use the mass code, 
	//          in the AtomTopology we usually have real masses
	at.setMass(d[0]-1);
        at.setCharge(d[1]);

        bs.addAtom(at);
    }

    // Constraints
    d_gin >> num;
    for (int j=0; j<num; j++){
      d_gin >> i[0] >> i[1] >> d[0];
      Constraint cons(--i[0],--i[1]);
      cons.setDist(d[0]);
      bs.addConstraint(cons);
    }
    if(! d_gin.check())
      throw InBuildingBlock::Exception("Building block file "+d_gin.name()+" is corrupted:\n MTBUILDBLSOLVENT block "+bs.solvName()+" is not OK!");
    d_bld.addBbSolvent(bs);
    d_gin >> s;
    if(s=="MTBUILDBLSOLVENT") cntn=1;
    else cntn=0;
  }
  // MTBUILDBLEND blocks
  cntn=0;
  if(!d_gin.eof()){
    if(s=="MTBUILDBLEND") cntn=1;
  }
  while(cntn){
    BbEnd be;
    int na, nrep;
    d_gin >> s;

    be.setResName(s);
    d_gin >> na >> nrep;
    be.setRep(nrep);
    for(int j = 0; j<na; j++){
	AtomTopology at;
	d_gin >> i[0];
        d_gin >> s;
        at.setName(s);
        d_gin >> i[0];
        at.setIac(--i[0]);
        d_gin >> d[0] >> d[1];
	// WARNING: in the building block we use the mass code, 
	//          in the AtomTopology we usually have real masses
	at.setMass(d[0]-1);
        at.setCharge(d[1]);
        d_gin >> i[0];
        at.setChargeGroup(i[0]);
        // The atoms that will replace atoms in a following residue do not
	// have exclusions specified. If nrep < 0, we are dealing with the
	// O-terminus.
        if(j<na - nrep){
	    Exclusion e;
	    d_gin >> i[0];
	    for(int k=0;k<i[0];k++){
		d_gin >> i[1];
		e.insert(--i[1]);
	    }
	    at.setExclusion(e);
	}
        be.addAtom(at);
    }

    // Bonds
    d_gin >> num;
    for (int j=0; j<num; j++){
      d_gin >> i[0] >> i[1] >> i[2];
      Bond bond(--i[0],--i[1]);
      bond.setType(--i[2]);
      be.addBond(bond);
    }

    // Angles
    d_gin >> num;
    for (int j=0; j<num; j++){
	d_gin >> i[0] >> i[1] >> i[2] >> i[3];
	Angle angle(--i[0],--i[1],--i[2]);
	angle.setType(--i[3]);
	be.addAngle(angle);
    }

    // Impropers
    d_gin >> num;
    for (int j=0; j<num; j++){
	d_gin >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
	Improper improper(--i[0],--i[1],--i[2],--i[3]);
	improper.setType(--i[4]);
	be.addImproper(improper);
    }

    // Dihedrals
    d_gin >> num;
    for (int j=0; j<num; j++){
	d_gin >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
	Dihedral dihedral(--i[0],--i[1],--i[2],--i[3]);
	dihedral.setType(--i[4]);
	be.addDihedral(dihedral);
    }


    if(! d_gin.check())
      throw InBuildingBlock::Exception("Building block file "+d_gin.name()+" is corrupted:\n MTBUILDBLEND block "+be.resName()+" is not OK!");
    d_bld.addBbEnd(be);
    d_gin >> s;
    if(s=="MTBUILDBLEND") cntn=1;
    else cntn=0; 
  }
}



















