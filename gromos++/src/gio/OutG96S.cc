// gio_OutG96S.cc

#include <cassert>
#include "OutG96S.h"
#include "../gromos/Exception.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"
#include "../gcore/Box.h"
#include <iostream>
#include <iomanip>

using gio::OutG96S;
using gio::OutG96S_i;
using namespace gcore;
using namespace std;

class OutG96S_i{
  friend class gio::OutG96S;
  ostream &d_os;
  int d_count, d_res_off, d_switch;
  OutG96S_i(ostream &os):
    d_os(os), d_count(0), d_res_off(1), d_switch()
  {d_switch = 0;}
  ~OutG96S_i(){}
  
  void select(const string &thing);
  void writeSingleM(const Molecule &mol);
  void writeSingleS(const Solvent &sol);
  void writeSingleM_vel(const Molecule &mol);
  void writeSingleS_vel(const Solvent &sol);
  void writeBox(const Box &box);  
  void writeTriclinicBox(const Box &box);
  void writeGenBox(const Box &box);
};

OutG96S::OutG96S(ostream &os):
  OutCoordinates(),
  d_this(new OutG96S_i(os)){}
OutG96S::OutG96S():
  OutCoordinates()
{d_this=0;}

OutG96S::~OutG96S(){
  if(d_this)delete d_this;
}

void OutG96S::writeTitle(const string &title){
  d_this->d_os << "TITLE\n" << title << "\nEND\n";
}

void OutG96S::writeTimestep(const int step, const double time)
{
  d_this->d_os.precision(9);
  d_this->d_os.setf(std::ios::fixed, std::ios::floatfield);

  d_this->d_os << "TIMESTEP\n"
               << std::setw(18)
               << step
               << std::setw(15)
               << time
               << "\nEND\n";
}

void OutG96S::select(const string &thing){
  if (thing == "ALL"){
      d_this->d_switch = 1;
  }
  else if (thing == "SOLVENT"){
    d_this->d_switch = 2;
  }
  else {
    d_this->d_switch = 0;
  }
}
  
void OutG96S::open(ostream &os){
  if(d_this){
    delete d_this;
  }
  d_this=new OutG96S_i(os);
}

void OutG96S::close(){
  if(d_this)delete d_this;
  d_this=0;
}

OutG96S &OutG96S::operator<<(const gcore::System &sys){
  // what do we have to write out
  bool writePos=false;
  bool writeVel=false;
  for(int m=0; m<sys.numMolecules(); m++){
    if(sys.mol(m).numPos()){
      writePos=true;
    }
    if(sys.mol(m).numVel()){
      writeVel=true;
    }
  }
  for(int s=0; s<sys.numSolvents(); s++){
    if(sys.sol(s).numPos()){
      writePos=true;
    }
    if(sys.sol(s).numVel()){
      writeVel=true;
    }
  }

  if(writePos){
    d_this->d_count=0;
    d_this->d_res_off=1;
    if(d_this->d_switch == 2){d_this->d_res_off=0;}
    d_this->d_os << "POSITION\n";
    if (d_this->d_switch <= 1)
      for(int i=0;i<sys.numMolecules();++i)
	d_this->writeSingleM(sys.mol(i));
    if (d_this->d_switch >= 1)
      for(int i=0;i<sys.numSolvents();++i)
	d_this->writeSingleS(sys.sol(i));
    
    d_this->d_os << "END\n";
  }

  if(writeVel){
    d_this->d_count=0;
    d_this->d_res_off=1;
    if(d_this->d_switch == 2){d_this->d_res_off=0;}
    d_this->d_os << "VELOCITY\n";
    if(d_this->d_switch <= 1){
      for(int i=0; i<sys.numMolecules();++i){
	d_this->writeSingleM_vel(sys.mol(i));
      }
    }
    if (d_this->d_switch >= 1){
      for(int i=0; i<sys.numSolvents(); ++i){
	d_this->writeSingleS_vel(sys.sol(i));
      }
    }
    d_this->d_os << "END\n";
  }
  switch(sys.box().boxformat()){
    case gcore::Box::box96:
      d_this->d_os << "BOX\n";
      d_this->writeBox(sys.box());
      d_this->d_os << "END\n";
      break;
    case gcore::Box::triclinicbox:
      d_this->d_os << "TRICLINICBOX\n";
      d_this->writeTriclinicBox(sys.box());
      d_this->d_os << "END\n";
      break;
    case gcore::Box::genbox:
      d_this->d_os << "GENBOX\n";
      d_this->writeGenBox(sys.box());
      d_this->d_os << "END\n";
      break;
    default:
      throw gromos::Exception("OutG96", "Don't know how to handle boxformat");
  }

  return *this;
}

void OutG96S_i::writeBox(const Box &box){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);

  d_os << setw(15) << box[0] 
       << setw(15) << box[1]
       << setw(15) << box[2] << endl;
}

void OutG96S_i::writeTriclinicBox(const Box &box){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);

  d_os << setw(8) << box.ntb() << endl;
  for(int i=0; i<3; ++i){
    d_os << setw(15) << box.K()[i] 
	 << setw(15) << box.L()[i]
	 << setw(15) << box.M()[i] << endl;
  }
  
}

void OutG96S_i::writeSingleM(const Molecule &mol){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(9);
  for (int i=0;i<mol.numAtoms();++i){
    ++d_count;
    int res=mol.topology().resNum(i);
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res + d_res_off;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' <<setw(6)<< mol.topology().resName(res).c_str()
	 << setw(6) << mol.topology().atom(i).name().c_str();
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(6) << d_count
	 << setw(15) << mol.pos(i)[0]
	 << setw(15) << mol.pos(i)[1]
	 << setw(15) << mol.pos(i)[2]<< endl;
  }
  d_res_off += mol.topology().numRes();
}

void OutG96S_i::writeSingleS(const Solvent &sol){
  int na=sol.topology().numAtoms();
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(9);
  for (int i=0;i<sol.numPos();++i){
    ++d_count;
    int res=i/na;
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res + d_res_off;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' <<setw(6)<< sol.topology().solvName().c_str()
	 << setw(6) << sol.topology().atom(i%na).name().c_str();
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(6) << d_count
	 << setw(15) << sol.pos(i)[0]
	 << setw(15) << sol.pos(i)[1]
	 << setw(15) << sol.pos(i)[2]<< endl;
  }
  d_res_off += sol.numPos()/na;
}

void OutG96S_i::writeSingleM_vel(const Molecule &mol){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(9);
  for (int i=0;i<mol.numVel();++i){
    ++d_count;
    int res=mol.topology().resNum(i);
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res + d_res_off;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' <<setw(6)<< mol.topology().resName(res).c_str()
         << setw(6) << mol.topology().atom(i).name().c_str();
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(6) << d_count
         << setw(15) << mol.vel(i)[0]
         << setw(15) << mol.vel(i)[1]
         << setw(15) << mol.vel(i)[2]<< endl;
  }
  d_res_off += mol.topology().numRes();
  
}

void OutG96S_i::writeSingleS_vel(const Solvent &sol){
  int na=sol.topology().numAtoms();
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(9);
  for (int i=0;i<sol.numVel();++i){
    ++d_count;
    int res=i/na;
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res + d_res_off;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' <<setw(6)<< sol.topology().solvName().c_str()
         << setw(6) << sol.topology().atom(i%na).name().c_str();
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(6) << d_count
         << setw(15) << sol.vel(i)[0]
         << setw(15) << sol.vel(i)[1]
         << setw(15) << sol.vel(i)[2]<< endl;
  }
  d_res_off += sol.numVel()/na;
}
void OutG96S_i::writeGenBox(const Box &box){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);
  const double k=box.K().abs();
  const double l=box.L().abs();
  const double m=box.M().abs();
  d_os << setw(8) << box.ntb() << endl;
  d_os << setw(15) << k
       << setw(15) << l
       << setw(15) << m << endl;
  d_os << setw(15) << acos(box.L().dot(box.M())/(l*m))*180/M_PI
       << setw(15) << acos(box.K().dot(box.M())/(k*m))*180/M_PI
       << setw(15) << acos(box.K().dot(box.L())/(k*l))*180/M_PI << endl;
  // construct a local x,y,z with x along k, y in the k,l plane and z in the direction of m

  Vec z = box.K().cross(box.L()).normalize();
  Vec x = box.K().normalize();
  Vec p,q;
  if(x[2]==0){
    p=x;
  }
  else{
    p=Vec(-z[1], z[0], 0);
    p=p.normalize();
  }
  q = -p.cross(z);
  
  double phi = acos (p.dot(x))*180/M_PI;
  double theta = asin(q[2])*180/M_PI;
  double psi = asin(p[1])*180/M_PI;

  d_os << setw(15) << phi 
       << setw(15) << theta
       << setw(15) << psi << endl;
}


