// gio_InBuildingBlock.cc

#include "InBuildingBlock.h"
#include "Ginstream.h"
#include "../gcore/BbSolute.h"
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
#include <sstream>

using namespace gcore;
using gio::InBuildingBlock_i;
using gio::InBuildingBlock;

// Implementation class
class InBuildingBlock_i: public gio::Ginstream
{
  friend class InBuildingBlock;
  gcore::BuildingBlock d_bld;
  /**
   * The init function reads in the whole file and stores
   * the information in the BuildingBlock
   */
  void init();
  /**
   * A function to read in MTBUILDBLSOLUTE blocks
   */
  void readSolute(std::vector<std::string> & buffer);
  void readEnd(std::vector<std::string> & buffer);
  void readSolvent(std::vector<std::string> &buffer);
  void readTopphyscon(std::vector<std::string> &buffer);
  void readLinkexclusions(std::vector<std::string> &buffer);
  
  InBuildingBlock_i (std::string &s): d_bld()
  {
    this->open(s);
    this->init();
  }
};

// Constructors

InBuildingBlock::InBuildingBlock(std::string name)
{
  d_this = new InBuildingBlock_i(name);
}

InBuildingBlock::~InBuildingBlock(){
  delete d_this;
}

const std::string InBuildingBlock::title()const{
  return d_this->title();
}

const BuildingBlock &InBuildingBlock::building()const{
  return d_this->d_bld;
}

void gio::InBuildingBlock_i::init()
{
  if(!stream())
    throw InBuildingBlock::Exception("Could not open Building block file "
				     +name());

  std::vector<std::string> buffer;

  while(!stream().eof()){
    getblock(buffer);
    if(buffer.size()){
      if     (buffer[0] == "MTBUILDBLSOLUTE")  readSolute(buffer);
      else if(buffer[0] == "MTBUILDBLSOLVENT") readSolvent(buffer);
      else if(buffer[0] == "MTBUILDBLEND")     readEnd(buffer);
      else if(buffer[0] == "TOPPHYSCON")       readTopphyscon(buffer);
      else if(buffer[0] == "LINKEXCLUSIONS")   readLinkexclusions(buffer);
    }
  }
}

void gio::InBuildingBlock_i::readTopphyscon(std::vector<std::string> &buffer)
{
  if(buffer.size() < 3) 
    throw InBuildingBlock::Exception("BuildingBlock file " + name() +
				     " is corrupted. Empty TOPPHYSCON block");
  std::string topphyscon;
  double d[2];
  
  gio::concatenate(buffer.begin()+1, buffer.end()-1, topphyscon);
  _lineStream.clear();
  _lineStream.str(topphyscon);
  _lineStream >> d[0] >> d[1];
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in TOPPHYSCON block:\n"+
				     topphyscon);
  d_bld.setFpepsi(d[0]);
  d_bld.setHbar(d[1]);
}

void gio::InBuildingBlock_i::readLinkexclusions(std::vector<std::string> &buffer)
{
  // The link exclusions are actually never used in the gromos++ way of 
  // making and editing topologies. This block is redundant.
  if(buffer.size() != 3) 
    throw InBuildingBlock::Exception("BuildingBlock file " + name() +
      " is corrupted. LINKEXCLUSIONS block should have only one line");
  int i;
  _lineStream.clear();
  _lineStream.str(buffer[1]);
  _lineStream >> i;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in LINKEXCLUSIONS block:\n"+
				     buffer[1]);
  d_bld.setLinkExclusions(i);
}

void gio::InBuildingBlock_i::readSolute(std::vector<std::string> &buffer)
{
  if(buffer.size()<3)
    throw InBuildingBlock::Exception("Building block file " + name()
		     + " is corrupted. Empty MTBUILDBLSOLUTE block.");
  // generic variables
  double d[4];
  int i[5], na, npe, num;
  std::string resname, s;
  
  std::string block;
  BbSolute bb;
  gio::concatenate(buffer.begin()+1, buffer.end()-1, block);
  
  _lineStream.clear();
  _lineStream.str(block);
  _lineStream >> resname;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block:\n"
			     +block+"Trying to read a residue name");
  bb.setResName(resname);

  _lineStream >> na >> npe;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block " 
			       + resname + ":\nTrying to read NMAT and NLIN");
  for(int j=0; j<npe; j++){
    Exclusion e;
    _lineStream >> i[0] >> i[1];
    if(i[0]!=j+1-npe)
      throw InBuildingBlock::Exception("Preceding exclusions in MTBUILDBLSOLUTE " + resname + " are not sequential");
    
    for(int k=0; k<i[1]; k++){
      _lineStream >> i[2];
      e.insert(--i[2]);
    }
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLUTE block " << resname 
	 << ".\nTrying to read preceding exclusion number " << j+1;
      throw InBuildingBlock::Exception(os.str());
    }
    bb.addPexcl(e);
  }
  
  for(int j = 0; j<na; j++){
    AtomTopology at;
    _lineStream >> i[0];
    if(i[0] != j+1)
      throw InBuildingBlock::Exception("Atom numbers in MTBUILDBLSOLUTE " + 
				      resname + " are not sequential");
    _lineStream >> s;
    at.setName(s);
    _lineStream >> i[0] >> d[0] >> d[1] >> i[1];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLUTE block " << resname
	 << ".\nTrying to read atom number " << j+1;
      throw InBuildingBlock::Exception(os.str());
    }
    at.setIac(--i[0]);
    // WARNING: in the building block we use the mass code, 
    //          in the AtomTopology we usually have real masses
    at.setMass(d[0]-1);
    at.setCharge(d[1]);
    at.setChargeGroup(i[1]);
    // The trailing atoms do not have exclusions specified.
    // these are specified as the preceding exclusions of the 
    // next residue
    if(j<na - npe){
      Exclusion e;
      _lineStream >> i[0];
      for(int k=0;k<i[0];k++){
	_lineStream >> i[1];
	e.insert(--i[1]);
      }
      if(_lineStream.fail()){
	std::ostringstream os;
	os << "Bad line in MTBUILDBLSOLUTE block " << resname 
	   << ".\nTrying to read exclusions of atom " << j+1;
	throw InBuildingBlock::Exception(os.str());
      }
      at.setExclusion(e);
    }
    bb.addAtom(at);
  }
  
  // Bonds
  _lineStream >> num;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block "
		+resname+".\nTrying to read number of bonds.");
  for (int j=0; j<num; j++){
    _lineStream >> i[0] >> i[1] >> i[2];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLUTE block " << resname
	 << ".\nTrying to read " << num << " bonds";
      throw InBuildingBlock::Exception(os.str());
    }
    Bond bond(--i[0],--i[1]);
    bond.setType(--i[2]);
    bb.addBond(bond);
  }
  // Angles
  _lineStream >> num;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block "
		+resname+".\nTrying to read number of angles.");
  for (int j=0; j<num; j++){
    _lineStream >> i[0] >> i[1] >> i[2] >> i[3];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLUTE block " << resname
	 << ".\nTrying to read " << num << " angles";
      throw InBuildingBlock::Exception(os.str());
    }
    Angle angle(--i[0],--i[1],--i[2]);
    angle.setType(--i[3]);
    bb.addAngle(angle);
  }
  // Impropers
  _lineStream >> num;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block "
				     +resname+".\nTrying to read number of impropers.");
  for (int j=0; j<num; j++){
    _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLUTE block " << resname
	 << ".\nTrying to read " << num << " impropers";
      throw InBuildingBlock::Exception(os.str());
    }
    Improper improper(--i[0],--i[1],--i[2],--i[3]);
    improper.setType(--i[4]);
    bb.addImproper(improper);
  }
  // Dihedrals
  _lineStream >> num;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block "
		    +resname+".\nTrying to read number of dihedrals.");
  for (int j=0; j<num; j++){
    _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLUTE block " << resname
	 << ".\nTrying to read " << num << " dihedrals";
      throw InBuildingBlock::Exception(os.str());
    }
    Dihedral dihedral(--i[0],--i[1],--i[2],--i[3]);
    dihedral.setType(--i[4]);
    bb.addDihedral(dihedral);
  }
  _lineStream >> s;
  if(!_lineStream.eof())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block "
		     + resname + ".\nTrailing data after dihedrals:" + s);
  d_bld.addBbSolute(bb);
}

void gio::InBuildingBlock_i::readSolvent(std::vector<std::string> &buffer)
{
  if(buffer.size()<3)
    throw InBuildingBlock::Exception("Building block file " + name()
		     + " is corrupted. Empty MTBUILDBLSOLUTE block.");
  // generic variables
  double d[4];
  int i[5], na, num;
  std::string resname, s;
  
  std::string block;
  gcore::SolventTopology st;
  gio::concatenate(buffer.begin()+1, buffer.end()-1, block);
  
  _lineStream.clear();
  _lineStream.str(block);
  _lineStream >> resname;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLVENT block:\n"
			     + block + "\nTrying to read the solvent name");
  st.setSolvName(resname);

  _lineStream >> na;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLVENT block " 
		     + resname + "\n:Trying to read the number of atoms");
  for(int j = 0; j<na; j++){
    AtomTopology at;
    _lineStream >> i[0];
    if(i[0] != j+1)
      throw InBuildingBlock::Exception("Atom numbers in MTBUILDBLSOLVENT " + 
				      resname + " are not sequential");
    _lineStream >> s;
    at.setName(s);
    _lineStream >> i[0] >> d[0] >> d[1];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLVENT block " << resname
	 << ".\nTrying to read atom number " << j+1;
      throw InBuildingBlock::Exception(os.str());
    }
    at.setIac(--i[0]);
    // WARNING: in the building block we use the mass code, 
    //          in the AtomTopology we usually have real masses
    at.setMass(d[0]-1);
    at.setCharge(d[1]);
    st.addAtom(at);
  }
  
  // Constraints
  _lineStream >> num;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLVENT block "
		+resname+".\nTrying to read number of constraints.");
  for (int j=0; j<num; j++){
    _lineStream >> i[0] >> i[1] >> d[0];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLVENT block " << resname
	 << ".\nTrying to read " << num << " constraints";
      throw InBuildingBlock::Exception(os.str());
    }
    Constraint cnstr(--i[0],--i[1]);
    cnstr.setDist(d[0]);
    st.addConstraint(cnstr);
  }
  _lineStream >> s;
  if(!_lineStream.eof())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block "
		     + resname + ".\nTrailing data after dihedrals:" + s);
  d_bld.addBbSolvent(st);
}

void gio::InBuildingBlock_i::readEnd(std::vector<std::string> &buffer)
{
  if(buffer.size()<3)
    throw InBuildingBlock::Exception("Building block file " + name()
		     + " is corrupted. Empty MTBUILDBLEND block.");
  // generic variables
  double d[4];
  int i[5], na, nrep, num;
  std::string resname, s;
  
  std::string block;
  gcore::BbSolute bb;
  gio::concatenate(buffer.begin()+1, buffer.end()-1, block);
  
  _lineStream.clear();
  _lineStream.str(block);
  _lineStream >> resname;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLEND block:\n"
			     + block + "\nTrying to read a residue name");
  bb.setResName(resname);

  _lineStream >> na >> nrep;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLEND block " 
       + resname + 
       ":\nTrying to read number of atoms and number of atoms to be replaced");
  bb.setRep(nrep);
  
  for(int j = 0; j<na; j++){
    AtomTopology at;
    _lineStream >> i[0];
    if(i[0] != j+1)
      throw InBuildingBlock::Exception("Atom numbers in MTBUILDBLSOLUTE " + 
				      resname + " are not sequential");
    _lineStream >> s;
    at.setName(s);
    _lineStream >> i[0] >> d[0] >> d[1] >> i[1];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLUTE block " << resname
	 << ".\nTrying to read atom number " << j+1;
      throw InBuildingBlock::Exception(os.str());
    }
    at.setIac(--i[0]);
    // WARNING: in the building block we use the mass code, 
    //          in the AtomTopology we usually have real masses
    at.setMass(d[0]-1);
    at.setCharge(d[1]);
    at.setChargeGroup(i[1]);
    // The atoms that will replace atoms in a following residue do not
    // have exclusions specified. If nrep < 0, we are dealing with the
    // O-terminus.not have exclusions specified.
    if(j<na - nrep){
      Exclusion e;
      _lineStream >> i[0];
      for(int k=0;k<i[0];k++){
	_lineStream >> i[1];
	e.insert(--i[1]);
      }
      if(_lineStream.fail()){
	std::ostringstream os;
	os << "Bad line in MTBUILDBLSOLUTE block " << resname 
	   << ".\nTrying to read exclusions of atom " << j+1;
	throw InBuildingBlock::Exception(os.str());
      }
      at.setExclusion(e);
    }
    bb.addAtom(at);
  }
  
  // Bonds
  _lineStream >> num;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block "
		+resname+".\nTrying to read number of bonds.");
  for (int j=0; j<num; j++){
    _lineStream >> i[0] >> i[1] >> i[2];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLUTE block " << resname
	 << ".\nTrying to read " << num << " bonds";
      throw InBuildingBlock::Exception(os.str());
    }
    Bond bond(--i[0],--i[1]);
    bond.setType(--i[2]);
    bb.addBond(bond);
  }
  // Angles
  _lineStream >> num;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block "
		+resname+".\nTrying to read number of angles.");
  for (int j=0; j<num; j++){
    _lineStream >> i[0] >> i[1] >> i[2] >> i[3];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLUTE block " << resname
	 << ".\nTrying to read " << num << " angles";
      throw InBuildingBlock::Exception(os.str());
    }
    Angle angle(--i[0],--i[1],--i[2]);
    angle.setType(--i[3]);
    bb.addAngle(angle);
  }
  // Impropers
  _lineStream >> num;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block "
				     +resname+".\nTrying to read number of impropers.");
  for (int j=0; j<num; j++){
    _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLUTE block " << resname
	 << ".\nTrying to read " << num << " impropers";
      throw InBuildingBlock::Exception(os.str());
    }
    Improper improper(--i[0],--i[1],--i[2],--i[3]);
    improper.setType(--i[4]);
    bb.addImproper(improper);
  }
  // Dihedrals
  _lineStream >> num;
  if(_lineStream.fail())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block "
		    +resname+".\nTrying to read number of dihedrals.");
  for (int j=0; j<num; j++){
    _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Bad line in MTBUILDBLSOLUTE block " << resname
	 << ".\nTrying to read " << num << " dihedrals";
      throw InBuildingBlock::Exception(os.str());
    }
    Dihedral dihedral(--i[0],--i[1],--i[2],--i[3]);
    dihedral.setType(--i[4]);
    bb.addDihedral(dihedral);
  }
  _lineStream >> s;
  if(!_lineStream.eof())
    throw InBuildingBlock::Exception("Bad line in MTBUILDBLSOLUTE block "
		     + resname + ".\nTrailing data after dihedrals:" + s);
  d_bld.addBbEnd(bb);
  
}



















