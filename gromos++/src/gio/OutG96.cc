// gio_OutG96.cc

#include "OutG96.h"
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

using gio::OutG96;
using gio::OutG96_i;
using namespace gcore;
using namespace std;

class OutG96_i{
  friend class OutG96;
  ostream &d_os;
  int d_count;
  int d_switch;
  OutG96_i(ostream &os):
    d_os(os), d_count(0)
  {d_switch = 0;}
  ~OutG96_i(){}

  void select(const string &thing);
  void writeTrajM(const Molecule &mol);
  void writeTrajS(const Solvent &sol);
  void writeBox(const Box &box);
};

OutG96::OutG96(ostream &os):
  OutCoordinates(),
  d_this(new OutG96_i(os)){}
OutG96::OutG96():
  OutCoordinates()
{d_this=0;}

OutG96::~OutG96(){
  if(d_this)delete d_this;
}

void OutG96::writeTitle(const string &title){
  d_this->d_os << "TITLE\n" << title << "\nEND\n";
}

void OutG96::select(const string &thing){
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

void OutG96::open(ostream &os){
  if(d_this){
    delete d_this;
  }
  d_this=new OutG96_i(os);
}

void OutG96::close(){
  if(d_this)delete d_this;
  d_this=0;
}

OutG96 &OutG96::operator<<(const gcore::System &sys){
  d_this->d_count=0;
  d_this->d_os << "POSITIONRED\n";
  if (d_this->d_switch <=1){
  for(int i=0;i<sys.numMolecules();++i){
    d_this->writeTrajM(sys.mol(i));}
  }
  if (d_this->d_switch >=1){
  for(int i=0;i<sys.numSolvents();++i){
    d_this->writeTrajS(sys.sol(i));}
  }
  d_this->d_os << "END\n";

  d_this->d_os << "BOX\n";
  d_this->writeBox(sys.box());
  d_this->d_os << "END\n";

  return *this;
}

void OutG96_i::writeTrajM(const Molecule &mol){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);
  for (int i=0;i<mol.numAtoms();++i){
    ++d_count;
    d_os << setw(15) << mol.pos(i)[0]
	 << setw(15) << mol.pos(i)[1]
	 << setw(15) << mol.pos(i)[2]<< endl;
    if(!(d_count%10))
      d_os << "#" << setw(10)<<d_count<<endl;
  } 
}
void OutG96_i::writeTrajS(const Solvent &sol){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);
  for (int i=0;i<sol.numCoords();++i){
    ++d_count;
    d_os << setw(15) << sol.pos(i)[0]
	 << setw(15) << sol.pos(i)[1]
	 << setw(15) << sol.pos(i)[2]<< endl;
    if(!(d_count%10))
      d_os << "#" << setw(10)<<d_count<<endl;
  } 
}


void OutG96_i::writeBox(const Box &box){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);

  d_os << setw(15) << box[0] 
       << setw(15) << box[1]
       << setw(15) << box[2] << endl;
}

