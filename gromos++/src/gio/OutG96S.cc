// gio_OutG96S.cc

#include "OutG96S.h"
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
  friend class OutG96S;
  ostream &d_os;
  int d_count, d_res_off, d_switch;
  OutG96S_i(ostream &os):
    d_os(os), d_count(0), d_res_off(1), d_switch()
  {d_switch = 0;}
  ~OutG96S_i(){}
  
  void select(const string &thing);
  void writeSingleM(const Molecule &mol);
  void writeSingleS(const Solvent &sol);
  void writeBox(const Box &box);
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

  d_this->d_os << "BOX\n";
  d_this->writeBox(sys.box());
  d_this->d_os << "END\n";

  return *this;
}

void OutG96S_i::writeBox(const Box &box){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);

  d_os << setw(15) << box[0] 
       << setw(15) << box[1]
       << setw(15) << box[2] << endl;
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
    d_os << setw(5) << d_count
	 << setw(15) << mol.pos(i)[0]
	 << setw(15) << mol.pos(i)[1]
	 << setw(15) << mol.pos(i)[2]<< endl;
  }
  d_res_off += mol.topology().resNum(mol.numAtoms()-1);
}

void OutG96S_i::writeSingleS(const Solvent &sol){
  
  int na=sol.topology().numAtoms();
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(9);
  for (int i=0;i<sol.numCoords();++i){
    ++d_count;
    int res=i/na + 1;
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res + d_res_off;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' <<setw(6)<< sol.topology().solvName().c_str()
	 << setw(6) << sol.topology().atom(i%na).name().c_str();
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << d_count
	 << setw(15) << sol.pos(i)[0]
	 << setw(15) << sol.pos(i)[1]
	 << setw(15) << sol.pos(i)[2]<< endl;
  }
  d_res_off += sol.numCoords()/na;
}


  
