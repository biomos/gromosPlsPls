// gio_OutPdb.cc

#include <cassert>
#include "OutPdb.h"
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

using gio::OutPdb;
using gio::OutPdb_i;
using namespace gcore;
using namespace std;

class OutPdb_i{
  friend class OutPdb;
  ostream &d_os;
  int d_count, d_resoff, d_switch;
  OutPdb_i(ostream &os):
    d_os(os), d_count(0), d_switch()
  {d_switch = 0;}
  ~OutPdb_i(){}

  void writeSingleM(const Molecule &mol);
  void writeSingleS(const Solvent &sol);
  void writeBox(const Box &box);
};

OutPdb::OutPdb(ostream &os):
  OutCoordinates(),
  d_this(new OutPdb_i(os)){}
OutPdb::OutPdb():
  OutCoordinates()
{d_this=0;}

OutPdb::~OutPdb(){
  if(d_this)delete d_this;
}

void OutPdb::writeTitle(const string &title){
  d_this->d_os << "TITLE " << title << "\n";
}

void OutPdb::select(const string &thing){
  if (thing == "ALL"){
    d_this->d_switch = 1;
  }
  else if (thing =="SOLVENT"){
    d_this->d_switch = 2;
  }
  else {
    d_this->d_switch = 0;
  }
}
  
void OutPdb::open(ostream &os){
  if(d_this){
    delete d_this;
  }
  d_this=new OutPdb_i(os);
}

void OutPdb::close(){
  if(d_this)delete d_this;
  d_this=0;
}

OutPdb &OutPdb::operator<<(const gcore::System &sys){
  d_this->d_os << "MODEL\n";
  d_this->d_count=0;
  d_this->d_resoff=1;
  if (d_this->d_switch <=1)
    for(int i=0;i<sys.numMolecules();++i)
      d_this->writeSingleM(sys.mol(i));
  if (d_this->d_switch >=1)
    for(int i=0;i<sys.numSolvents();++i)
      d_this->writeSingleS(sys.sol(i));
  d_this->d_os << "TER\n";
  d_this->d_os << "ENDMDL\n";

  return *this;
}

void OutPdb_i::writeBox(const Box &box){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);

  d_os << setw(15) << box[0] 
       << setw(15) << box[1]
       << setw(15) << box[2] << endl;
}

void OutPdb_i::writeSingleM(const Molecule &mol){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);
  for (int i=0; i<mol.numAtoms(); ++i){
    ++d_count;
    int res=mol.topology().resNum(i);
    d_os << "ATOM";
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(7) << d_count;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << "  " <<setw(4) << mol.topology().atom(i).name().substr(0,3).c_str();
    d_os << setw(4) << mol.topology().resName(res).substr(0,4).c_str();
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res+d_resoff << "    "
	 << setw(8) << mol.pos(i)[0]*10
	 << setw(8) << mol.pos(i)[1]*10
	 << setw(8) << mol.pos(i)[2]*10
	 << "  1.00  0.00" << endl;
  }
  d_resoff+=mol.topology().numRes();
  
}
void OutPdb_i::writeSingleS(const Solvent &sol){
  int na=sol.topology().numAtoms();
  
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);
  for (int i=0; i<sol.numPos(); ++i){
    ++d_count;
    int res=i/na;
    
    d_os << "ATOM";
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(7) << d_count;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << "  " <<setw(4) << sol.topology().atom(i%na).name().substr(0,3).c_str();
    d_os << setw(4) << sol.topology().solvName().substr(0,4).c_str();
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res+d_resoff << "    "
	 << setw(8) << sol.pos(i)[0]*10
	 << setw(8) << sol.pos(i)[1]*10
	 << setw(8) << sol.pos(i)[2]*10
	 << "  1.00  0.00" << endl;
  }
  d_resoff+=sol.numPos()/na;
  
}
