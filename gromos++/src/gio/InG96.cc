// gio_InG96.cc

#include "InG96.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Vec.h"
#include "Ginstream.h"
#include "../gcore/Solvent.h"
#include "../gcore/Box.h"
#include <map>

enum blocktype {timestep,positionred,position,velocity,box};

typedef map<string,blocktype>::value_type BT;
// Define class variable BT (block_types)
const BT blocktypes[] = {BT("TIMESTEP",timestep),
			 BT("POSITIONRED", positionred),
			 BT("POSITION", position),
			 BT("VELOCITY",velocity),
			 BT("BOX", box)};

static map<string,blocktype> BLOCKTYPE(blocktypes,blocktypes+5);

using gio::InG96;
using gio::InG96_i;
using namespace gcore;

class InG96_i{
  friend class InG96;
  Ginstream d_gin;
  string d_name;
  string d_current;
  int d_switch;
  InG96_i(const string &name):
    d_gin(name), d_name(name), d_current(), d_switch(){d_gin >> d_current;
    d_switch = 0;}
  ~InG96_i(){}
  // method
  void readMolecule(Molecule &);
  void readSolvent(Solvent &, int);
  void readBox(Box &);
};

// Constructors
InG96::InG96(const string &name): d_this(new InG96_i(name)){}
InG96::InG96(){d_this=0;}

InG96::~InG96(){
  if(d_this)delete d_this;
}

void InG96::open(const string &name){
  if(d_this)delete d_this;
  d_this=new InG96_i(name);
}

void InG96::close(){
  if(d_this)delete d_this;
  d_this=0;
}

void InG96::select(const string &thing){
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
  return d_this->d_gin.eof();
}

const string &InG96::title()const{
  return d_this->d_gin.title();
}

void InG96_i::readMolecule(Molecule &m){
  int na=m.numAtoms();

  if(d_current=="POSITION"){
    string dummy;
    for (int i=0;i<na;i++){
      for(int j=0;j<4;j++) d_gin >> dummy;
      d_gin >> m.pos(i)[0] >> m.pos(i)[1] >> m.pos(i)[2];
    } 
    return;
  }
  
  else{
      for (int i=0;i<na;i++)
	d_gin >> m.pos(i)[0] >> m.pos(i)[1] >> m.pos(i)[2];
      return;
  }
  throw InG96::Exception("POSITION or POSITIONRED block is corrupted in "+d_name);
      
}

void InG96_i::readSolvent(Solvent &s, int dum){
  s.setnumCoords(0);
  Vec v; 

  if(d_current=="POSITION"){
    string dummy;
    // get rid of solute
    if(d_switch == 2){
   for (int i=0;i<dum;++i){
     d_gin.getline(dummy,100);}}
    d_gin >> dummy;
    while(dummy!="END"){
      for(int j=0;j<3;j++) d_gin >> dummy;
      d_gin >> v[0] >> v[1] >> v[2];
      d_gin >> dummy;
      s.addCoord(v);
    }
    return;
  }
 
  else{
   string dummy;
     // get rid of solute
     if(d_switch == 2){
         for (int i=0;i<dum;++i){
	   d_gin.getline(dummy,100);}}  
       d_gin >> dummy;
    while(dummy!="END"){
      v[0]=atof(dummy.c_str());
      d_gin >> v[1] >> v[2];
      d_gin >> dummy;
      s.addCoord(v);  
    }   
   return;
  }
  throw InG96::Exception("POSITION or POSITIONRED block is corrupted in "+d_name);

  }

void InG96_i::readBox(Box &box){
  d_gin >> box[0] >> box[1] >> box[2];
}

InG96 &InG96::operator>>(System &sys){

  if(!d_this->d_gin)
    throw Exception("File "+name()+" is corrupted.");
  
  string first;

  first=d_this->d_current;
  int nm=sys.numMolecules();
  int dum=0;
   for (int i=0;i<nm;++i){
  dum += sys.mol(i).numAtoms();
  } 

  
  do{
    switch(BLOCKTYPE[d_this->d_current]){
    case timestep:
      // not yet implemented
      while(! d_this->d_gin.check());
      break;
    case positionred:
      if (d_this->d_switch <=1){ 
        for(int i=0;i<nm;i++){
          d_this->readMolecule(sys.mol(i));}}
      if (d_this->d_switch >=1){
        d_this->readSolvent(sys.sol(0), dum);}
      else{
        while(! d_this->d_gin.check());}
      break;
    case velocity:
      // not yet implemented
      while(! d_this->d_gin.check());
      break;
    case position:
       if (d_this->d_switch <=1){ 
    for(int i=0;i<nm;i++){
      d_this->readMolecule(sys.mol(i));}}
       if (d_this->d_switch >=1){
	 d_this->readSolvent(sys.sol(0), dum);}
            else{
        while(! d_this->d_gin.check());}
      break;
    case box:
      d_this->readBox(sys.box());
      if(! d_this->d_gin.check())
	throw Exception("BOX block is not OK in "+name()+".");
      break;
    default:
      throw
	Exception("Block "+d_this->d_current+" is unknown in a coordinates file");
      break;
    }
    d_this->d_gin>>d_this->d_current;
  } while(d_this->d_current!=first&&!d_this->d_gin.eof());
  return *this;
  
}

const string &InG96::name()const{
  return d_this->d_name;
}
