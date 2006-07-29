// gio_InG96.cc

#include <cassert>
#include "InG96.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Vec.h"
#include "../gmath/Matrix.h"
#include "Ginstream.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Box.h"
#include <map>

enum blocktype {timestep,
		positionred,position,
		velocityred,velocity,
		box,triclinicbox, genbox, gmxbox};

typedef std::map<std::string,blocktype>::value_type BT;
// Define class variable BT (block_types)
const BT blocktypes[] = {BT("TIMESTEP",timestep),
			 BT("POSITIONRED", positionred),
			 BT("POSITION", position),
			 BT("VELOCITYRED", velocityred),
			 BT("VELOCITY",velocity),
			 BT("BOX", box),
                         BT("TRICLINICBOX", triclinicbox),
			 BT("GENBOX", genbox),
                         BT("GMXBOX", gmxbox)};

static std::map<std::string,blocktype> BLOCKTYPE(blocktypes,blocktypes+9);

using gio::InG96;
using gio::InG96_i;
using namespace gcore;

class InG96_i: public gio::Ginstream
{
  friend class gio::InG96;

  std::string d_current;
  int d_switch;

  int d_skip;
  int d_stride;
  
  InG96_i(const std::string &name, int skip, int stride)
    : d_current(),
      d_switch(),
      d_skip(skip),
      d_stride(stride)
  {
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
  void readGmxbox(gcore::System &sys);
  void readBox(gcore::System &sys);
  void readGenbox(gcore::System &sys);
};

// Constructors
InG96::InG96(const std::string &name, int skip, int stride)
  : d_skip(skip),
    d_stride(stride),
    d_stride_eof(false)
{
  d_this=new InG96_i(name, skip, stride);
}

InG96::InG96(int skip, int stride)
  : d_this(NULL),
    d_skip(skip),
    d_stride(stride),
    d_stride_eof(false)
{
}

InG96::~InG96(){
  if(d_this)delete d_this;
}

void InG96::open(const std::string &name){
  if(d_this){
    // recover skip and stride
    d_skip = d_this->d_skip;
    d_stride = d_this->d_stride;
    
    delete d_this;
  }
  
  d_this=new InG96_i(name, d_skip, d_stride);
}

void InG96::close(){
  if(d_this){
    d_this->close();

    d_skip = d_this->d_skip;
    d_stride = d_this->d_stride;
    
    delete d_this;
    d_this=NULL;
  }
}

void InG96::select(const std::string &thing){
  if (!d_this){
    throw InG96::Exception("select must be called after open!");
  }
  
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
  if (!d_this){
    throw InG96::Exception("eof, but no open file");
  }
  
  return d_this->stream().eof();
}

std::string InG96::title()const{
  if (!d_this){
    throw InG96::Exception("no open file");
  }
  return d_this->title();
}
			 
void InG96_i::readTimestep(gcore::System &sys)
{
  std::vector<std::string> buffer;
  getblock(buffer);
  // not implemented;
}

void InG96_i::readPosition(gcore::System &sys)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;
  getblock(buffer);
  if(buffer[buffer.size()-1].find("END")!=0)
    throw InG96::Exception("Coordinate file " + name() +
			   " is corrupted. No END in POSITION"
			   " block. Got\n"
			   + buffer[buffer.size()-1]);

  it=buffer.begin();
  std::string dummy;
  int begin=0;
  if(d_current=="POSITION") begin=24;
  int na=0;
  for(int m=0; m<sys.numMolecules(); m++){
    if(!sys.mol(m).numPos()){
      sys.mol(m).initPos();
    }
    na+=sys.mol(m).numAtoms();
  }
  
  // solute?
  if(d_switch<2){
    for(int m=0; m<sys.numMolecules(); m++){
      for(int a=0; a<sys.mol(m).numAtoms(); a++, ++it){
	if(unsigned (begin)>=(*it).size()){
	  std::ostringstream os;
	  os << "Coordinate file " << name() << " corrupted.\n"
	     << "Failed to read coordinates from line \n"
	     << *it
	     << "\nfrom POSITION or POSITIONRED block";
	  throw InG96::Exception(os.str());
	}
	
	_lineStream.clear();
	_lineStream.str((*it).substr(begin, (*it).size()));
	_lineStream >> sys.mol(m).pos(a)[0]
		    >> sys.mol(m).pos(a)[1]
		    >> sys.mol(m).pos(a)[2];
	
	
	if(_lineStream.fail()){
	  std::ostringstream os;
	  os << "Coordinate file " << name() << " corrupted.\n"
	     << "Failed to read " << na << " solute coordinates"
	     << " from POSITION or POSITIONRED block";
	  throw InG96::Exception(os.str());
	}
      }
      
    }
  } else it+=na;
  
  // Solvent?
  if(d_switch > 0){
    sys.sol(0).setNumPos(0);
    gmath::Vec v;
    
    for(; it!=buffer.end()-1; ++it){
      if(unsigned(begin)>=(*it).size()){
	std::ostringstream os;
	os << "Coordinate file " << name() << " corrupted.\n"
	   << "Failed to read coordinates from line \n"
	   << *it
	   << "\nfrom POSITION or POSITIONRED block";
	throw InG96::Exception(os.str());
      }
      _lineStream.clear();
      _lineStream.str((*it).substr(begin, (*it).size()));
      _lineStream >> v[0] >> v[1] >> v[2];
      sys.sol(0).addPos(v);
    
      if(_lineStream.fail()){
	std::ostringstream os;
	os << "Coordinate file " << name() << " corrupted.\n"
	   << "Failed while reading solvent coordinates"
	   << " from POSITION or POSITIONRED block";
	throw InG96::Exception(os.str());
      }
    }
    
    if(sys.sol(0).numPos() % sys.sol(0).topology().numAtoms() != 0){
      std::ostringstream os;
      os << "Coordinate file " << name() << " corrupted.\n"
	 << "Atom count mismatch while reading solvent coordinates.\n"
	 << "Read " << sys.sol(0).numPos() << " coordinates for solvent "
	 << "with " << sys.sol(0).topology().numAtoms() << " atoms per molecule\n";
      throw InG96::Exception(os.str());
    }
  }
}

void InG96_i::readVelocity(gcore::System &sys)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;
  getblock(buffer);
  if(buffer[buffer.size()-1].find("END")!=0)
    throw InG96::Exception("Coordinate file " + name() +
			   " is corrupted. No END in VELOCITY"
			   " block. Got\n"
			   + buffer[buffer.size()-1]);

  it=buffer.begin();
  std::string dummy;
  int begin=0;
  if(d_current=="VELOCITY") begin=24;
  int na=0;
  for(int m=0; m<sys.numMolecules(); m++){
    if(!sys.mol(m).numVel()) {
      sys.mol(m).initVel();
    }
    na+=sys.mol(m).numVel();
  }
  // solute?
  if(d_switch<2){
    for(int m=0; m<sys.numMolecules(); m++){
      for(int a=0; a<sys.mol(m).numVel(); a++, ++it){ 
	if(unsigned (begin) >=(*it).size()){
	  std::ostringstream os;
	  os << "Coordinate file " << name() << " corrupted.\n"
	     << "Failed to read velocities from line \n"
	     << *it
	     << "\nfrom VELOCITY or VELOCITYRED block";
	  throw InG96::Exception(os.str());
	}
        _lineStream.clear();
        _lineStream.str((*it).substr(begin, (*it).size()));
        _lineStream >> sys.mol(m).vel(a)[0]
                    >> sys.mol(m).vel(a)[1]
                    >> sys.mol(m).vel(a)[2];
    
	if(_lineStream.fail()){
	  std::ostringstream os;
	  os << "Coordinate file " << name() << " corrupted.\n"
	     << "Failed to read " << na << " solute velocity coordinates"
	     << " from VELOCITY or VELOCITYRED block";
	  throw InG96::Exception(os.str());
	}
      }
    }
  } else it+=na;
  
  // Solvent?
  if(d_switch > 0){
    sys.sol(0).setNumVel(0);
    gmath::Vec v;
    
    for(; it!=buffer.end()-1; ++it){
      if(unsigned (begin) >=(*it).size()){
	std::ostringstream os;
	os << "Coordinate file " << name() << " corrupted.\n"
	   << "Failed to read velocities from line \n"
	   << *it
	   << "\nfrom VELOCITTY or VELOCITYRED block";
	throw InG96::Exception(os.str());
      }
      _lineStream.clear();
      _lineStream.str((*it).substr(begin, (*it).size()));
      _lineStream >> v[0] >> v[1] >> v[2];
      sys.sol(0).addVel(v);

      if(_lineStream.fail()){
	std::ostringstream os;
	os << "Coordinate file " << name() << " corrupted.\n"
	   << "Failed while reading solvent velocity coordinates"
	   << " from VELOCITY or VELOCITYRED block";
	throw InG96::Exception(os.str());
      }
    }
    if(sys.sol(0).numVel() % sys.sol(0).topology().numAtoms() != 0){
      std::ostringstream os;
      os << "Coordinate file " << name() << " corrupted.\n"
	 << "Atom count mismatch while reading solvent velocities.\n"
	 << "Read " << sys.sol(0).numVel() << " coordinates for solvent "
	 << "with " << sys.sol(0).topology().numAtoms() << " atoms per molecule\n";
      throw InG96::Exception(os.str());
    }
  }
}

void InG96_i::readBox(gcore::System &sys){
  // std::cerr << "readbox" << std::endl;
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;
  getblock(buffer);
  if(buffer[buffer.size()-1].find("END")!=0)
    throw InG96::Exception("Coordinate file " + name() +
			   " is corrupted. No END in BOX"
			   " block. Got\n"
			   + buffer[buffer.size()-1]);
  
  it=buffer.begin();
  
  _lineStream.clear();
  _lineStream.str(*it);
  _lineStream >> sys.box()[0] >> sys.box()[1] >> sys.box()[2];
  if(_lineStream.fail())
    throw InG96::Exception("Bad line in BOX block:\n" + *it + 
			   "\nTrying to read three doubles");  
}

void InG96_i::readGmxbox(System &sys)
{   
  std::vector<std::string> buffer;                            
  std::vector<std::string>::iterator it;                      
  
  getblock(buffer);                                           
  if(buffer[buffer.size()-1].find("END")!=0)                  
    throw InG96::Exception("Coordinate file " + name() +      
                           " is corrupted. No END in GMXBOX"
                           " block. Got\n"                    
                           + buffer[buffer.size()-1]);        
  
  std::string s;
  gio::concatenate(buffer.begin(), buffer.end()-1, s);        
                                                              
  _lineStream.clear();                                        
  _lineStream.str(s);                                         
  sys.box().setNtb(gcore::Box::boxshape_enum(gcore::Box::triclinic));           

  gmath::Vec k,l,m;                                           
  _lineStream >> k[0] >> l[1] >> m[2]                         
              >> k[1] >> k[2] >> l[0]                         
              >> l[2] >> m[0] >> m[1];                        
  sys.box().K()=k;                                            
  sys.box().L()=l;                                            
  sys.box().M()=m;                                            
                                                              
  sys.box().update_triclinic();                               
                                                              
  if(_lineStream.fail())                                      
    throw InG96::Exception("Bad line in GMXBOX block:\n" + *it +      
                           "\nTrying to read nine doubles");  
                                                              
}
void InG96_i::readTriclinicbox(System &sys)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;
  
  getblock(buffer);
  if(buffer[buffer.size()-1].find("END")!=0)
    throw InG96::Exception("Coordinate file " + name() +
			   " is corrupted. No END in TRICLINICBOX"
			   " block. Got\n"
			   + buffer[buffer.size()-1]);

  std::string s;
  gio::concatenate(buffer.begin(), buffer.end()-1, s);
  
  _lineStream.clear();
  _lineStream.str(s);
  int ntb;
  _lineStream >> ntb;

  // GromosXX truncated octahedral box:  3
  // Gromos++ truncated octahedral box: -1
  // yes, i know!
  //if (ntb == -1) ntb = 3;

  sys.box().setNtb(gcore::Box::boxshape_enum(ntb));

  gmath::Vec k,l,m;
  _lineStream >> k[0] >> l[0] >> m[0]
	      >> k[1] >> l[1] >> m[1]
	      >> k[2] >> l[2] >> m[2];
  sys.box().K()=k;
  sys.box().L()=l;
  sys.box().M()=m;

  sys.box().update_triclinic();
  
  if(_lineStream.fail())
    throw InG96::Exception("Bad line in TRICLINICBOX block:\n" + *it + 
			   "\nTrying to read nine doubles");

}

void InG96_i::readGenbox(System &sys)
{
  std::vector<std::string> buffer;
  
  getblock(buffer);
  if(buffer[buffer.size()-1].find("END")!=0)
    throw InG96::Exception("Coordinate file " + name() +
			   " is corrupted. No END in GENBOX"
			   " block. Got\n"
			   + buffer[buffer.size()-1]);

  std::string s;
  gio::concatenate(buffer.begin(), buffer.end()-1, s);
  
  _lineStream.clear();
  _lineStream.str(s);
  int ntb;
  _lineStream >> ntb;
  sys.box().setNtb(gcore::Box::boxshape_enum(ntb));
  sys.box().boxformat() = gcore::Box::genbox;
  double a,b,c,alpha,beta,gamma,phi,theta,psi;
  _lineStream >> a >> b >> c
	      >> alpha >> beta >> gamma 
	      >> phi >> theta >> psi;
  
  if(_lineStream.fail())
    throw InG96::Exception("Bad line in GENBOX block:\n" + s + 
			   "\nTrying to read nine doubles");

  if(ntb==gcore::Box::vacuum)
    return;
  if(ntb==gcore::Box::rectangular || ntb==gcore::Box::truncoct){
    if(alpha!=90 || beta!=90 || gamma!=90)
      throw InG96::Exception("For rectangular and truncated octahedral boxes, alpha, beta"
			     " and gamma should be 90 degrees");
    if(phi!=0 || theta != 0 || psi !=0)
      throw InG96::Exception("For rectangular and truncated octahedral boxes, phi, theta"
			     " and phi should be 0 degrees");
    sys.box().K()=Vec(a, 0.0,0.0);
    sys.box().L()=Vec(0.0,b ,0.0);
    sys.box().M()=Vec(0.0,0.0,c );
    return;
  }
  alpha *=M_PI/180.0;
  beta  *=M_PI/180.0;
  gamma *=M_PI/180.0;
  phi   *=-M_PI/180.0;
  theta *=-M_PI/180.0;
  psi   *=-M_PI/180.0;
  
  const double cosphi   =cos(phi);
  const double sinphi   =sin(phi);
  const double costheta =cos(theta);
  const double sintheta =sin(theta);
  const double cospsi   =cos(psi);
  const double sinpsi   =sin(psi);

  sys.box().K() = Vec(1,0,0);
  sys.box().L() = Vec(cos(gamma), sin(gamma),0);
  const double cosbeta=cos(beta);
  const double cosalpha=cos(alpha);
  // const double dotp=sys.box().L()[1];

  // double mx=cosbeta;
  // double my=cosalpha*dotp;
  
  sys.box().M() = Vec(cosbeta, 1,
		      sqrt((cosbeta*sys.box().L()[0] + sys.box().L()[1]) * 
			   (cosbeta*sys.box().L()[0] + sys.box().L()[1]) / (cosalpha * cosalpha)
			   - cosbeta*cosbeta - 1) );
  
			   

  //sys.box().M() = Vec(cosbeta,
  //		      -(cosalpha - sys.box().M()[0] * sys.box().L()[0] / sys.box().L()[1]),
  //		      sqrt(1 - cosbeta*cosbeta - (cosalpha - sys.box().M()[0] * sys.box().L()[0] / sys.box().L()[1]) * (cosalpha - sys.box().M()[0] * sys.box().L()[0] / sys.box().L()[1])));

  sys.box().M() = sys.box().M().normalize();
  sys.box().K() *=a;
  sys.box().L() *=b;
  sys.box().M() *=c;
  
  Vec x(cospsi*cosphi - costheta*sinphi*sinpsi,
	-sinpsi*cosphi - costheta*sinphi*cospsi,
	sintheta*sinphi);
  Vec y(cospsi*sinphi + costheta*cosphi*sinpsi,
	-sinpsi*sinphi + costheta*cosphi*cospsi,
	-sintheta*cosphi);
  Vec z(sinpsi*sintheta,
	cospsi*sintheta,
	costheta);
  gmath::Matrix mat(x,y,z);
  sys.box().K() = mat*sys.box().K();
  sys.box().L() = mat*sys.box().L();
  sys.box().M() = mat*sys.box().M();
  std::cout << gmath::v2s(x) << std::endl;
  std::cout << gmath::v2s(y) << std::endl;
  std::cout << gmath::v2s(z) << std::endl;
  /**
   * This seems to work with the wrong conventions
   sys.box().K() = x;
  sys.box().L() = x*cos(gamma) + y * sin(gamma);
  const double cosbeta=cos(beta);
  const double cosalpha=cos(alpha);
  const double dotp=sys.box().L().dot(y);
  
  sys.box().M() = cosbeta * x + cosalpha * dotp*y 
    + sqrt(1 -cosbeta * cosbeta - cosalpha * cosalpha * dotp * dotp) * z;

  sys.box().K() *= a;
  sys.box().L() *= b;
  sys.box().M() *= c;
  */
}

InG96 &InG96::operator>>(System &sys){

  if (!d_this){
    throw InG96::Exception("read in frame, but no open file");
  }
  
  if(!d_this->stream())
    throw Exception("File "+name()+" is corrupted.");
  
  const std::string first =d_this->d_current;
  // std::cerr << first << std::endl;
  std::vector<std::string> buffer;
  bool readpos = false;
  bool readvel = false;
  bool readbox = false;

  // skip frames
  // std::cerr << "operator<< : skip=" << d_this->d_skip << std::endl;
  
  for( ; d_this->d_skip > 0; --d_this->d_skip){
    
    do{
      // std::cerr << "skipping block " << d_this->d_current << std::endl;
      
      d_this->skipblock();
      d_this->getline(d_this->d_current);
      
    } while (d_this->d_current != first &&
	     (!d_this->stream().eof()));
    
    if (d_this->stream().eof()){
      // std::cerr << "skip eof: " << d_this->d_skip << std::endl;

      --d_this->d_skip;
      d_stride_eof = true;
      return *this;
    }
    
  }

  // only stride if not skip because of eof during last stride
  if (d_stride_eof == false){
    int i=1;
    for( ; i < d_this->d_stride; ++i){
      
      do{
	
	d_this->skipblock();
	d_this->getline(d_this->d_current);
	
      } while (d_this->d_current != first &&
	       (!d_this->stream().eof()));
      
      if (d_this->stream().eof()){
	// safe remaining strides in skip for next file
	
	std::cerr << "stride eof: " << d_this->d_stride
		  << "\ti: " << i << std::endl;
	
	d_this->d_skip = d_this->d_stride - i - 1;
	d_stride_eof = true;
	return *this;
      }
      
    }
  }
  
  // skipping and striding worked...
  d_stride_eof = false;

  do{
    switch(BLOCKTYPE[d_this->d_current]){
      case timestep:
	// not yet implemented
	d_this->readTimestep(sys);
	break;
      case positionred:
	d_this->readPosition(sys);
	readpos = true;
	break;
      case position:
	d_this->readPosition(sys);
	readpos = true;
	break;
      case velocityred:
	d_this->readVelocity(sys);
        readvel = true;
	break;
      case velocity:
	d_this->readVelocity(sys);
        readvel = true;
	break;
      case box:
	d_this->readBox(sys);
        readbox = true;
	break;
      case triclinicbox:
	d_this->readTriclinicbox(sys);
        readbox = true;
	break;
      case gmxbox:
        d_this->readGmxbox(sys);
        readbox = true;
        break;
      case genbox:
	d_this->readGenbox(sys);
        readbox = true;
	break;
      default:
	throw
	  Exception("Block "+d_this->d_current+
		    " is unknown in a coordinate file");
	break;
    }
    d_this->getline(d_this->d_current);
  } while(d_this->d_current!=first&&!d_this->stream().eof());

  sys.hasPos = readpos;
  sys.hasVel = readvel;
  sys.hasBox = readbox;
  return *this;
}

std::string InG96::name()const{
  return d_this->name();
}

int InG96::skip()const
{
  if (d_this)
    return d_this->d_skip;
  else return d_skip;
}

int InG96::stride()const
{
  if (d_this)
    return d_this->d_stride;
  else return d_stride;
}

void InG96::skip(int skip)
{
  if (d_this)
    d_this->d_skip = skip;
  else
    d_skip = skip;
}

void InG96::stride(int stride)
{
  if (d_this)
    d_this->d_stride = stride;
  else
    d_stride = stride;
}

bool InG96::stride_eof()const
{
  return d_stride_eof;
}
