// gio_OutTopology.cc

#include "OutTopology.h"
#include "../gcore/BondType.h"
#include "../gcore/Bond.h"
#include "../gcore/AngleType.h"
#include "../gcore/Angle.h"
#include "../gcore/DihedralType.h"
#include "../gcore/Dihedral.h"
#include "../gcore/ImproperType.h"
#include "../gcore/Improper.h"
#include "../gcore/LJType.h"
#include "../gcore/AtomPair.h"
#include "../gcore/Exclusion.h"
#include "../gcore/Constraint.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/System.h"
#include "../gcore/GromosForceField.h"
#include <iostream>
#include <iomanip>

using gio::OutTopology;
using namespace gcore;
using namespace std;

OutTopology::OutTopology(ostream &os):d_os(os){
}

OutTopology::~OutTopology(){}

void OutTopology::setTitle(const string &title){
  d_title=title;
}

//OutTopology &OutTopology::operator<<(const gcore::Simulation &sim){
void OutTopology::write(const gcore::System &sys, const gcore::GromosForceField &gff){
  
  // Title block
  d_os << "TITLE\n" << d_title << "\nEND\n";

  // TOPPHYSCON block
  d_os.precision(10);
  
  d_os << "TOPPHYSCON\n" << gff.fpepsi() << endl << gff.hbar()<<"\nEND\n";
  
  // TOPVERSION block
  d_os << "TOPVERSION\n1.7\nEND\n";

  // ATOMTYPENAME block
  d_os << "ATOMTYPENAME\n";
  int num=gff.numAtomTypeNames();
  d_os << num << endl;
  for(int i=0;i<num;++i){
    if(i>0 &&!(i%10))d_os << "# " << i << endl;
    d_os << gff.atomTypeName(i) << endl;
  }
  d_os << "END\n";
  
  // RESNAME block
  d_os << "RESNAME\n";
  num=0;
  for(int i=0;i<sys.numMolecules();++i)
    num+=sys.mol(i).topology().numRes();

  d_os << num << endl;
  for(int i=0, count=0;i<sys.numMolecules();++i)
    for(int j=0;j<sys.mol(i).topology().numRes();++j,++count){
      if(count>0 &&!(count%10))d_os << "# " << count << endl;
      d_os << sys.mol(i).topology().resName(j) << endl;
  }
  d_os << "END\n";
  
  // SOLUTEATOM block
  d_os << "SOLUTEATOM" << endl;
  num=0;
  for(int i=0;i<sys.numMolecules();++i)
    num+=sys.mol(i).numAtoms();

  d_os << setw(5) << num << endl;

  for(int i=0, offatom=1, offres=1;i<sys.numMolecules();++i){
    for(int j=0;j<sys.mol(i).numAtoms();++j){
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(5) << offatom+j << ' '
	   << setw(5) << sys.mol(i).topology().resNum(j)+offres << ' '
	   << setw(4) << sys.mol(i).topology().atom(j).name()
	   << setw(4) << sys.mol(i).topology().atom(j).iac()+1
	   << setw(11) << sys.mol(i).topology().atom(j).mass()
	   << setw(11) << sys.mol(i).topology().atom(j).charge() 
	   << setw(3) << sys.mol(i).topology().atom(j).chargeGroup();
      // Exclusions
      d_os << setw(3) << sys.mol(i).topology().atom(j).exclusion().size();
      for(int k=0;k<sys.mol(i).topology().atom(j).exclusion().size();++k){
	if(!((k+1)%6))d_os << endl << "                                             ";
	d_os << setw(6) << sys.mol(i).topology().atom(j).exclusion().atom(k)+offatom;
      }
      d_os << endl;
      d_os << "                                             " 
	   << sys.mol(i).topology().atom(j).exclusion14().size();
      for(int k=0;k<sys.mol(i).topology().atom(j).exclusion14().size();++k)
	d_os << setw(6) << sys.mol(i).topology().atom(j).exclusion14().atom(k)+offatom;
      d_os << endl;
    }
    offres+=sys.mol(i).topology().numRes();
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // BONDTYPE block

  d_os << "BONDTYPE\n";
  num=gff.numBondTypes();

  d_os << num << endl;
  for (int i=0;i<num;++i){
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.bondType(i).fc()
	 << setw(12) << gff.bondType(i).b0() << endl;
  }
  d_os << "END\n";

  // BONDH block
  d_os << "BONDH\n";
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH())
	++num;
    }
  }
  d_os << num << endl;
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH()) 
	d_os << setw(6) << bit()[0] +offatom
	     << setw(6) << bit()[1]+offatom
	     << setw(5) << bit().type()+1 << endl;
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  d_os << "BOND\n";
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH()) 
	++num;
    }
  }
  d_os << num << endl;
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH()) 
	d_os << setw(6) << bit()[0] +offatom
	     << setw(6) << bit()[1]+offatom
	     << setw(5) << bit().type()+1 << endl;
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // BONDANGLETYPE block

  d_os << "BONDANGLETYPE\n";
  num=gff.numAngleTypes();

  d_os << num << endl;
  for (int i=0;i<num;++i){
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.angleType(i).fc()
	 << setw(12) << gff.angleType(i).t0() << endl;
  }
  d_os << "END\n";

  // BONDANGLEH & BONDANGLE block
  d_os << "BONDANGLEH\n";
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    AngleIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH())
	++num;
    }
  }
  d_os << num << endl;
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    AngleIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH())
	d_os << setw(6) << bit()[0] +offatom
	     << setw(6) << bit()[1]+offatom
	     << setw(6) << bit()[2]+ offatom
	     << setw(5) << bit().type()+1 << endl;
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  d_os << "BONDANGLE\n";
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    AngleIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH())
	  ++num;
    }
  }
  d_os << num << endl;
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    AngleIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH())   
	d_os << setw(6) << bit()[0] +offatom
	     << setw(6) << bit()[1]+offatom
	     << setw(6) << bit()[2] +offatom
	     << setw(5) << bit().type()+1 << endl;
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";
      
  // IMPDIHEDRALTYPE block

  d_os << "IMPDIHEDRALTYPE\n";
  num=gff.numImproperTypes();

  d_os << num << endl;
  for (int i=0;i<num;++i){
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.improperType(i).fc()
	 << setw(12) << gff.improperType(i).q0() << endl;
  }
  d_os << "END\n";

  // IMPDIHEDRALH & IMPDIHEDRAL block
  d_os << "IMPDIHEDRALH\n";
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    ImproperIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH() ||
	 sys.mol(i).topology().atom(bit()[3]).isH())
	++num;
    }
  }
  d_os << num << endl;
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    ImproperIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH() ||
	 sys.mol(i).topology().atom(bit()[3]).isH())
	d_os << setw(6) << bit()[0] +offatom
	     << setw(6) << bit()[1]+offatom
	     << setw(6) << bit()[2]+ offatom
	     << setw(6) << bit()[3]+ offatom
	     << setw(5) << bit().type()+1 << endl;
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  d_os << "IMPDIHEDRAL\n";
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    ImproperIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH() &&
	 !sys.mol(i).topology().atom(bit()[3]).isH())
	++num;
    }
  }
  d_os << num << endl;
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    ImproperIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH() &&
	 !sys.mol(i).topology().atom(bit()[3]).isH())
	d_os << setw(6) << bit()[0] +offatom
	     << setw(6) << bit()[1]+offatom
	     << setw(6) << bit()[2] +offatom
	     << setw(6) << bit()[3] + offatom
	     << setw(5) << bit().type()+1 << endl;
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // DIHEDRALTYPE block

  d_os << "DIHEDRALTYPE\n";
  num=gff.numDihedralTypes();

  d_os << num << endl;
  for (int i=0;i<num;++i){
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os << setw(10) << gff.dihedralType(i).fc()
	 << setw(10) << gff.dihedralType(i).pd() 
	 << setw(4)<< gff.dihedralType(i).np() << endl;
  }
  d_os << "END\n";

  // DIHEDRALH & DIHEDRAL block
  d_os << "DIHEDRALH\n";
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    DihedralIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH() ||
	 sys.mol(i).topology().atom(bit()[3]).isH())	
	++num;
    }
  }
  d_os << num << endl;
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    DihedralIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH() ||
	 sys.mol(i).topology().atom(bit()[3]).isH())	
	d_os << setw(6) << bit()[0] +offatom
	     << setw(6) << bit()[1]+offatom
	     << setw(6) << bit()[2]+ offatom
	     << setw(6) << bit()[3]+ offatom
	     << setw(5) << bit().type()+1 << endl;
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  d_os << "DIHEDRAL\n";
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    DihedralIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH() &&
	 !sys.mol(i).topology().atom(bit()[3]).isH())
	++num;
    }
  }
  d_os << num << endl;
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    DihedralIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH() &&
	 !sys.mol(i).topology().atom(bit()[3]).isH())
	d_os << setw(6) << bit()[0] +offatom
	     << setw(6) << bit()[1]+offatom
	     << setw(6) << bit()[2] +offatom
	     << setw(6) << bit()[3] + offatom
	     << setw(5) << bit().type()+1 << endl;
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";
  
  // LJPARAMETERS block
  d_os << "LJPARAMETERS\n";
  num=gff.numLJTypes();
  d_os << num << endl;
  for (int i=0;i<gff.numAtomTypeNames();++i)
    for(int j=0;j<=i;++j){
      d_os.precision(6);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os.setf(ios::scientific, ios::floatfield);
      LJType lj(gff.ljType(AtomPair(i,j)));
      d_os << setw(5) << j+1 
	   << setw(5) << i+1
	   << setw(14) << lj.c12()
	   << setw(14) << lj.c6()
	   << setw(14) << lj.cs12()
	   << setw(14) << lj.cs6() << endl;
    }

  d_os << "END\n";

  //SOLVENTATOM block
  d_os << "SOLVENTATOM\n";
  d_os << sys.sol(0).topology().numAtoms() << endl;

  for(int j=0;j<sys.sol(0).topology().numAtoms();++j){
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os << setw(5) << j+1
	 << setw(5) << sys.sol(0).topology().atom(j).name()
	 << setw(4) << sys.sol(0).topology().atom(j).iac()+1
	 << setw(11) << sys.sol(0).topology().atom(j).mass()
	 << setw(11) << sys.sol(0).topology().atom(j).charge();
    d_os << endl;
  }
  d_os << "END\n";

  //SOLVENTCONSTR bock
  d_os << "SOLVENTCONSTR\n";
  num=0;
  ConstraintIterator dit(sys.sol(0).topology());
  for(;dit;++dit) ++num;
  d_os << num << endl;
  
  ConstraintIterator cit(sys.sol(0).topology());

  for(;cit;++cit){
    d_os.precision(7);
    d_os.setf(ios::fixed, ios::floatfield);

    d_os << setw(5) << cit()[0] + 1
         << setw(5) << cit()[1] + 1
         << setw(15) << cit().dist() << endl;
  }
  d_os << "END\n";

  //  return *this;

}


