// gio_InG96.cc

#include <cassert>
#include "InG96.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Vec.h"
#include "Ginstream.h"
#include "../gcore/Solvent.h"
#include "../gcore/Box.h"
#include <map>

enum blocktype {timestep,
		positionred,position,
		velocityred,velocity,
		box,triclinicbox};

typedef std::map<std::string,blocktype>::value_type BT;
// Define class variable BT (block_types)
const BT blocktypes[] = {BT("TIMESTEP",timestep),
			 BT("POSITIONRED", positionred),
			 BT("POSITION", position),
			 BT("VELOCITYRED", velocityred),
			 BT("VELOCITY",velocity),
			 BT("BOX", box),
                         BT("TRICLINICBOX", triclinicbox)};

static std::map<std::string,blocktype> BLOCKTYPE(blocktypes,blocktypes+7);

using gio::InG96;
using gio::InG96_i;
using namespace gcore;

class InG96_i: public gio::Ginstream
{
  friend class InG96;
  std::string d_current;
  int d_switch;
  InG96_i(const std::string &name):
    d_current(), d_switch(){
    open(name);
    getline(d_current);
    d_switch = 0;
  }
  ~InG96_i(){}
  // method
  void readTimestep(gcore::System &sys);
  void readPosition(gcore::System &sys);
  void readVelocity(gcore::System &sys);
  void readTriclinicbox(gcore::System &sys);
  void readBox(gcore::System &sys);
};

// Constructors
InG96::InG96(const std::string &name){
  d_this=new InG96_i(name);
}

InG96::InG96(){
  d_this=0;
}

InG96::~InG96(){
  if(d_this)delete d_this;
}

void InG96::open(const std::string &name){
  if(d_this)delete d_this;
  d_this=new InG96_i(name);
}

void InG96::close(){
  d_this->close();
  if(d_this)delete d_this;
  d_this=0;
}

void InG96::select(const std::string &thing){
  if (thing == "ALL"){
    d_this->d_switch = 1;
  }
  else if (thing =="SOLVENT"){
    d_this->d_switch = 2;
  }  
  else if (thing =="SOLUTE"){
    d_this->d_switch = 0;
  }
  else {
    d_this->d_switch = 0;
  }
}

bool InG96::eof()const{
  return d_this->stream().eof();
}

const std::string InG96::title()const{
  // std::cerr << "readtitle" << std::endl;
  return d_this->title();
}
			 
void InG96_i::readTimestep(gcore::System &sys)
{
  // std::cerr << "readtime" << std::endl;
  std::vector<std::string> buffer;
  getblock(buffer);
  // not implemented;
}

void InG96_i::readPosition(gcore::System &sys)
{
  // std::cerr << "readpos" << std::endl;
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;
  getblock(buffer);
  it=buffer.begin();
  std::string dummy;
  int begin=0;
  if(d_current=="POSITION") begin=24;
  int na=0;
  for(int m=0; m<sys.numMolecules(); m++)
    na+=sys.mol(0).numAtoms();
  
  // solute?
  if(d_switch<2){
    for(int m=0; m<sys.numMolecules(); m++){
      for(int a=0; a<sys.mol(m).numAtoms(); a++, ++it){
	_lineStream.clear();
	_lineStream.str((*it).substr(begin, (*it).size()));
	_lineStream >> sys.mol(m).pos(a)[0]
		    >> sys.mol(m).pos(a)[1]
		    >> sys.mol(m).pos(a)[2];
      }
    
      if(_lineStream.fail()){
	std::ostringstream os;
	os << "Coordinate file " << name() << " corrupted.\n"
	   << "Failed to read " << na << " solute coordinates"
	   << " from POSITION or POSITIONRED block";
	throw InG96::Exception(os.str());
      }
      
    }
  } else it+=na;
  
  // Solvent?
  if(d_switch > 0){
    sys.sol(0).setnumCoords(0);
    gmath::Vec v;
    
    for(; it!=buffer.end()-1; ++it){
      _lineStream.clear();
      _lineStream.str((*it).substr(begin, (*it).size()));
      _lineStream >> v[0] >> v[1] >> v[2];
      sys.sol(0).addCoord(v);
    }
    if(_lineStream.fail()){
      std::ostringstream os;
      os << "Coordinate file " << name() << " corrupted.\n"
	 << "Failed while reading solvent coordinates"
	 << " from POSITION or POSITIONRED block";
      throw InG96::Exception(os.str());
    }
  }
}

void InG96_i::readVelocity(gcore::System &sys)
{
  std::vector<std::string> buffer;
  getblock(buffer);
  // not implemented
}

void InG96_i::readBox(gcore::System &sys){
  // std::cerr << "readbox" << std::endl;
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;
  getblock(buffer);
  it=buffer.begin();
  
  _lineStream.clear();
  _lineStream.str(*it);
  _lineStream >> sys.box()[0] >> sys.box()[1] >> sys.box()[2];
  if(_lineStream.fail())
    throw InG96::Exception("Bad line in BOX block:\n" + *it + 
			   "\nTrying to read three doubles");
}

void InG96_i::readTriclinicbox(System &sys)
{
  std::vector<std::string> buffer;
  getblock(buffer);
  //not implemented
}

InG96 &InG96::operator>>(System &sys){

  if(!d_this->stream())
    throw Exception("File "+name()+" is corrupted.");
  
  const std::string first =d_this->d_current;
  // std::cerr << first << std::endl;
  std::vector<std::string> buffer;
  
  do{
    switch(BLOCKTYPE[d_this->d_current]){
      case timestep:
	// not yet implemented
	d_this->readTimestep(sys);
	break;
      case positionred:
	d_this->readPosition(sys);
	break;
      case position:
	d_this->readPosition(sys);
	break;
      case velocityred:
	d_this->readVelocity(sys);
	break;
      case velocity:
	d_this->readVelocity(sys);
	break;
      case box:
	d_this->readBox(sys);
	break;
      case triclinicbox:
	d_this->readTriclinicbox(sys);
	break;
      default:
	throw
	  Exception("Block "+d_this->d_current+
		    " is unknown in a coordinate file");
	break;
    }
    d_this->getline(d_this->d_current);
  } while(d_this->d_current!=first&&!d_this->stream().eof());
  return *this;
}

const std::string InG96::name()const{
  return d_this->name();
}
