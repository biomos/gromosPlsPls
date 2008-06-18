// gio_OutTopology.cc

#include <cassert>
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
#include "../args/Arguments.h"
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


  if(args::Arguments::outG96){
    // TOPPHYSCON block
    d_os.precision(10);

    d_os << "TOPPHYSCON\n"
         << "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)\n"
         << gff.fpepsi()
         << "\n# HBAR: Planck's constant HBAR = H/(2* PI)\n"
         << gff.hbar()
         <<"\nEND\n";
  }
  else{
    // PHYSICALCONSTANTS block
    d_os.precision(10);
  
    d_os << "PHYSICALCONSTANTS\n"
         << "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)\n"
         << gff.fpepsi() 
         << "\n# HBAR: Planck's constant HBAR = H/(2* PI)\n"
         << gff.hbar()
         << "\n# BOLTZ: Boltzmann's constant kB\n"
         << gff.boltz()
         <<"\nEND\n";
  }
  
  // TOPVERSION block
  if(args::Arguments::outG96){
    d_os << "TOPVERSION\n1.7\nEND\n";
  }
  else{
    d_os << "TOPVERSION\n2.0\nEND\n";
  }

  // ATOMTYPENAME block
  d_os << "ATOMTYPENAME\n"
       << "# NRATT: number of van der Waals atom types\n";
  int num=gff.numAtomTypeNames();
  d_os << num << "\n";
  d_os << "# TYPE: atom type names\n";
  for(int i=0;i<num;++i){
    if(i>0 &&!(i%10))d_os << "# " << i << "\n";
    d_os << gff.atomTypeName(i) << "\n";
  }
  d_os << "END\n";
  
  // RESNAME block
  d_os << "RESNAME\n"
       << "# NRAA2: number of residues in a solute molecule\n";
  num=0;
  for(int i=0;i<sys.numMolecules();++i)
    num+=sys.mol(i).topology().numRes();

  d_os << num << "\n"
       << "# AANM: residue names\n";
  for(int i=0, count=0;i<sys.numMolecules();++i)
    for(int j=0;j<sys.mol(i).topology().numRes();++j,++count){
      if(count>0 &&!(count%10))d_os << "# " << count << "\n";
      d_os << sys.mol(i).topology().resName(j) << "\n";
  }
  d_os << "END\n";
  
  // SOLUTEATOM block
  d_os << "SOLUTEATOM\n"
       << "#   NRP: number of solute atoms\n";
  num=0;
  for(int i=0;i<sys.numMolecules();++i)
    num+=sys.mol(i).numAtoms();

  d_os << setw(5) << num << "\n";
  d_os << "#  ATNM: atom number\n"
       << "#  MRES: residue number\n"
       << "#  PANM: atom name of solute atom\n"
       << "#   IAC: integer (van der Waals) atom type code\n"
       << "#  MASS: mass of solute atom\n"
       << "#    CG: charge of solute atom\n"
       << "#   CGC: charge group code (0 or 1)\n"
       << "#   INE: number of excluded atoms\n"
       << "# INE14: number of 1-4 interactions\n"
       << "# ATNM MRES PANM IAC     MASS       CG  CGC INE\n"
       << "#                                           INE14\n";

  for(int i=0, offatom=1, offres=1;i<sys.numMolecules();++i){
    for(int j=0;j<sys.mol(i).numAtoms();++j){
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(6) << offatom+j << ' '
	   << setw(4) << sys.mol(i).topology().resNum(j)+offres << ' '
	   << setw(4) << sys.mol(i).topology().atom(j).name()
	   << setw(4) << sys.mol(i).topology().atom(j).iac()+1
	   << setw(9) << sys.mol(i).topology().atom(j).mass()
	   << setw(9) << sys.mol(i).topology().atom(j).charge() 
	   << setw(3) << sys.mol(i).topology().atom(j).chargeGroup();
      // Exclusions
      d_os << setw(3) << sys.mol(i).topology().atom(j).exclusion().size();
      for(int k=0;k<sys.mol(i).topology().atom(j).exclusion().size();++k){
	if(k%6==0 && k!=0)
	  d_os << "\n"
	       << "                                            ";
	d_os << setw(6) << sys.mol(i).topology().atom(j).exclusion().atom(k)+offatom;
      }
      d_os << "\n"
	   << "                                           " 
	   << sys.mol(i).topology().atom(j).exclusion14().size();
      for(int k=0;k<sys.mol(i).topology().atom(j).exclusion14().size();++k){
	if(k%6==0 && k!=0)
	  d_os << "\n"
	       << "                                            ";
	d_os << setw(6) << sys.mol(i).topology().atom(j).exclusion14().atom(k)+offatom;
      }
      
      d_os << "\n";
    }
    offres+=sys.mol(i).topology().numRes();
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // BONDTYPE block
  d_os << "BONDTYPE\n"
       << "#  NBTY: number of covalent bond types\n";
  num=gff.numBondTypes();

  d_os << num << "\n"
       << "#  CB: force constant\n"
       << "#  B0: bond length at minimum energy\n"
       << "#         CB          B0\n";

  for (int i=0;i<num;++i){
    if(i>0 &&!(i%10)) d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.bondType(i).fc()
     << setw(12) << gff.bondType(i).b0() << "\n";
  }
  d_os << "END\n";

  if(!(args::Arguments::outG96)){
    // BONDSTRETCHTYPE block

    d_os << "BONDSTRETCHTYPE\n"
         << "#  NBTY: number of covalent bond types\n";
    num=gff.numBondTypes();
  
    d_os << num << "\n"
         << "#  CB:  quartic force constant\n"
         << "#  CHB: harmonic force constant\n"
         << "#  B0:  bond length at minimum energy\n"
         << "#         CB         CHB         B0\n";

    for (int i=0;i<num;++i){
      if(i>0 &&!(i%10)) d_os << "# " << i << "\n";
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os.setf(ios::scientific, ios::floatfield);
      d_os << setw(16) << gff.bondType(i).fc()
           << setw(16) << gff.bondType(i).hfc()
	       << setw(16) << gff.bondType(i).b0() << "\n";
    }
    d_os << "END\n";
  }
  
  // BONDH block
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH())
	++num;
    }
  }
  d_os << "BONDH\n"
       << "#  NBONH: number of bonds involving H atoms in solute\n"
       << num << "\n"
       << "#  IBH, JBH: atom sequence numbers of atoms forming a bond\n"
       << "#  ICBH: bond type code\n"
       << "#   IBH    JBH ICBH\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";
  
  d_os << "BOND\n"
       << "#  NBON: number of bonds NOT involving H atoms in solute\n";
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH()) 
	++num;
    }
  }
  d_os << num << "\n"
       << "#  IB, JB: atom sequence numbers of atoms forming a bond\n"
       << "#  ICB: bond type code\n"
       << "#    IB     JB  ICB\n";
  
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH()) {
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      } 
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // BONDANGLETYPE block
  num=gff.numAngleTypes();
  d_os << "BONDANGLETYPE\n"
       << "#  NTTY: number of bond angle types\n"
       << num << "\n"
       << "#  CT: force constant\n"
       << "#  T0: bond angle at minimum energy in degrees\n"
       << "#         CT          T0\n";

  for (int i=0;i<num;++i){
    if(i>0 &&!(i%10))d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.angleType(i).fc()
     << setw(12) << gff.angleType(i).t0() << "\n";
  }
  d_os << "END\n";

  if(!(args::Arguments::outG96)){
    // BONDANGLEBENDTYPE block
    num=gff.numAngleTypes();
    d_os << "BONDANGLEBENDTYPE\n"
         << "#  NTTY: number of bond angle types\n"
         << num << "\n"
         << "#  CT:  force constant (based on potential\n"
         << "#       harmonic in the angle cosine)\n"
         << "#  CHT: force constant (based on potential\n"
         << "#       harmonic in the angle)\n"
         << "#  T0:  bond angle at minimum energy in degrees\n"
         << "#         CT         CHT          T0\n";

    for (int i=0;i<num;++i){
      if(i>0 &&!(i%10))d_os << "# " << i << "\n";
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os.setf(ios::scientific, ios::floatfield);
      d_os << setw(16) << gff.angleType(i).fc()
           << setw(16) << gff.angleType(i).afc()
	       << setw(16) << gff.angleType(i).t0() << "\n";
    }
    d_os << "END\n";
  }

  // BONDANGLEH & BONDANGLE block
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
  d_os << "BONDANGLEH\n"
       << "#  NTHEH: number of bond angles involving H atoms in solute\n"
       << num << "\n"
       << "#  ITH, JTH, KTH: atom sequence numbers\n"
       << "#    of atoms forming a bond angle in solute\n"
       << "#  ICTH: bond angle type code\n"
       << "#   ITH    JTH    KTH ICTH\n";
  
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    AngleIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH()){
	
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2]+ offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

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
  d_os << "BONDANGLE\n"
       << "#  NTHE: number of bond angles NOT\n"
       << "#        involving H atoms in solute\n"
       << num << "\n"
       << "#  IT, JT, KT: atom sequence numbers of atoms\n"
       << "#     forming a bond angle\n"
       << "#  ICT: bond angle type code\n"
       << "#    IT     JT     KT  ICT\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    AngleIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2] +offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";
      
  // IMPDIHEDRALTYPE block
  num=gff.numImproperTypes();
  d_os << "IMPDIHEDRALTYPE\n"
       << "#  NQTY: number of improper dihedrals\n"
       << num << "\n"
       << "#  CQ: force constant of improper dihedral per degrees square\n"
       << "#  Q0: improper dihedral angle at minimum energy in degrees\n"
       << "#         CQ          Q0\n";

  for (int i=0;i<num;++i){
    if(i>0 &&!(i%10))d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.improperType(i).fc()
	 << setw(12) << gff.improperType(i).q0() << "\n";
  }
  d_os << "END\n";

  // IMPDIHEDRALH & IMPDIHEDRAL block
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
  d_os << "IMPDIHEDRALH\n"
       << "#  NQHIH: number of improper dihedrals\n"
       << "#         involving H atoms in the solute\n"
       << num << "\n"
       << "#  IQH,JQH,KQH,LQH: atom sequence numbers\n"
       << "#     of atoms forming an improper dihedral\n"
       << "#  ICQH: improper dihedral type code\n"
       << "#   IQH    JQH    KQH    LQH ICQH\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    ImproperIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH() ||
	 sys.mol(i).topology().atom(bit()[3]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2]+ offatom
	     << setw(7) << bit()[3]+ offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

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
  d_os << "IMPDIHEDRAL\n"
       << "#  NQHI: number of improper dihedrals NOT\n"
       << "#    involving H atoms in solute\n"
       << num << "\n"
       << "#  IQ,JQ,KQ,LQ: atom sequence numbers of atoms\n"
       << "#    forming an improper dihedral\n"
       << "#  ICQ: improper dihedral type code\n"
       << "#    IQ     JQ     KQ     LQ  ICQ\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    ImproperIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH() &&
	 !sys.mol(i).topology().atom(bit()[3]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2] +offatom
	     << setw(7) << bit()[3] + offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // DIHEDRALTYPE block
  num=gff.numDihedralTypes();

  d_os << "DIHEDRALTYPE\n"
       << "#  NPTY: number of dihedral types\n"
       << num << "\n"
       << "#  CP: force constant\n"
       << "#  PD: cosine of the phase shift\n"
       << "#  NP: multiplicity\n"
       << "#       CP        PD  NP\n";

  for (int i=0;i<num;++i){
    if(i>0 &&!(i%10))d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os << setw(10) << gff.dihedralType(i).fc()
     << setw(10) << gff.dihedralType(i).pd()
     << setw(4)<< gff.dihedralType(i).np() << "\n";
  }
  d_os << "END\n";

  if(!(args::Arguments::outG96)){
    // TORSDIHEDRALTYPE block
    num=gff.numDihedralTypes();

    d_os << "TORSDIHEDRALTYPE\n"
         << "#  NPTY: number of dihedral types\n"
         << num << "\n"
         << "#  CP: force constant\n"
         << "#  PD: phase-shift angle\n"
         << "#  NP: multiplicity\n"
         << "#       CP        PD  NP\n";

    for (int i=0;i<num;++i){
      if(i>0 &&!(i%10)) d_os << "# " << i << "\n";
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(10) << gff.dihedralType(i).fc()
	       << setw(10) << gff.dihedralType(i).pdl() 
	       << setw(4)<< gff.dihedralType(i).np() << "\n";
    }
    d_os << "END\n";
  }

  // DIHEDRALH & DIHEDRAL block
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
  d_os << "DIHEDRALH\n"
       << "#  NPHIH: number of dihedrals involving H atoms in solute\n"
       << num << "\n"
       << "#  IPH, JPH, KPH, LPH: atom sequence numbers\n"
       << "#    of atoms forming a dihedral\n"
       << "#  ICPH: dihedral type code\n"
       << "#   IPH    JPH    KPH    LPH ICPH\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    DihedralIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH() ||
	 sys.mol(i).topology().atom(bit()[3]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2]+ offatom
	     << setw(7) << bit()[3]+ offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

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
  d_os << "DIHEDRAL\n"
       << "#  NPHI: number of dihedrals NOT involving H atoms in solute\n"
       << num << "\n"
       << "#  IP, JP, KP, LP: atom sequence numbers\n"
       << "#     of atoms forming a dihedral\n"
       << "#  ICP: dihedral type code\n"
       << "#    IP     JP     KP     LP  ICP\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    DihedralIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH() &&
	 !sys.mol(i).topology().atom(bit()[3]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2] +offatom
	     << setw(7) << bit()[3] + offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";
  
  // LJPARAMETERS block
  num=gff.numLJTypes();
  d_os << "LJPARAMETERS\n"
       << "#  NRATT2: number of LJ interaction types = NRATT*(NRATT+1)/2\n"
       << num << "\n"
       << "#  IAC,JAC: integer (van der Waals) atom type code\n"
       << "#  C12: r**(-12) term in nonbonded interactions\n"
       << "#   C6: r**(-6) term in nonbonded interactions\n"
       << "# CS12: r**(-12) term in 1-4 nonbonded interactions\n"
       << "#  CS6: r**(-6) term in 1-4 nonbonded interactions\n"
       << "# IAC  JAC           C12            C6          CS12           CS6\n";

  for (int i=0;i<gff.numAtomTypeNames();++i){
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
	   << setw(14) << lj.cs6() << "\n";
    }
    d_os << "#\n";
  }
  d_os << "END\n";


  if(!(args::Arguments::outG96)){
    // SOLUTEMOLECULES block
    // Still implemented for default only
    d_os << "SOLUTEMOLECULES\n"
         << "# NSPM: number of separate molecules in solute block\n"
         << "# NSP[1...NSPM]: atom sequence number of last atom\n"
         << "#                of the successive submolecules\n"
         << "#      NSPM  NSP[1...NSPM]\n";
  
    d_os << setw(10) << sys.numMolecules() << "\n";
    int nspmin = 0;
    for (int i=0; i<sys.numMolecules(); ++i){
      d_os << setw(6) << sys.mol(i).numAtoms()+nspmin;
      nspmin+=sys.mol(i).numAtoms();
      if((i+1)%10==0) d_os << "\n";
    }
 
    if(sys.numMolecules()%10!=0) d_os << "\n"; 
    d_os << "END\n";

    // TEMPERATUREGROUPS block
    // Still implemented for default only
    d_os << "TEMPERATUREGROUPS\n"
         << "# NSTM: number of temperature atom groups\n"
         << "# NST[1...NSTM]: atom sequence number of last atom\n"
         << "#                of the successive temperature atom groups\n"
         << "#      NSTM  NST[1...NSTM]\n";

    d_os << setw(10) << sys.numTemperatureGroups() << "\n";
    for(int i=0;i<sys.numTemperatureGroups();++i){
      d_os << setw(6)  << sys.temperatureGroup(i);
      if((i+1)%10==0) d_os << "\n";
    }

    if(sys.numTemperatureGroups()%10!=0) d_os << "\n";
    d_os << "END\n";
 
    // PRESSUREGROUPS block
    // Still implemented for default only
    d_os << "PRESSUREGROUPS\n"
         << "# NSVM: number of pressure atom groups\n"
         << "# NSV[1...NSVM]: atom sequence number of last atom\n"
         << "#                of the successive pressure atom groups\n"
         << "#      NSVM  NSV[1...NSVM]\n";

    d_os << setw(10) << sys.numPressureGroups() << "\n";
    for(int i=0;i<sys.numPressureGroups();++i){
      d_os << setw(6)  << sys.pressureGroup(i); 
      if((i+1)%10==0) d_os << "\n";
    }

    if(sys.numPressureGroups()%10!=0) d_os << "\n";
    d_os << "END\n";
  }

  //SOLVENTATOM block
  d_os << "SOLVENTATOM\n"
       << "#  NRAM: number of atoms per solvent molecule\n"
       << sys.sol(0).topology().numAtoms() << "\n"
       << "#     I: solvent atom sequence number\n"
       << "#  IACS: integer (van der Waals) atom type code\n"
       << "#  ANMS: atom name of solvent atom\n"
       << "#  MASS: mass of solvent atom\n"
       << "#   CGS: charge of solvent atom\n"
       << "#  I  ANMS IACS      MASS        CGS\n";

  for(int j=0;j<sys.sol(0).topology().numAtoms();++j){
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os << setw(4) << j+1 << " "
	 << setw(5) << sys.sol(0).topology().atom(j).name()
	 << setw(4) << sys.sol(0).topology().atom(j).iac()+1
	 << setw(11) << sys.sol(0).topology().atom(j).mass()
	 << setw(11) << sys.sol(0).topology().atom(j).charge();
    d_os << "\n";
  }
  d_os << "END\n";

  //SOLVENTCONSTR bock
  num=0;
  ConstraintIterator dit(sys.sol(0).topology());
  for(;dit;++dit) ++num;
  d_os << "SOLVENTCONSTR\n"
       << "#  NCONS: number of constraints\n"
       << num << "\n"
       << "#  ICONS, JCONS: atom sequence numbers forming constraint\n"
       << "#   CONS constraint length\n"
       << "#ICONS JCONS         CONS\n";

  ConstraintIterator cit(sys.sol(0).topology());

  for(;cit;++cit){
    d_os.precision(7);
    d_os.setf(ios::fixed, ios::floatfield);

    d_os << setw(5) << cit()[0] + 1
         << setw(5) << cit()[1] + 1
         << setw(15) << cit().dist() << "\n";
  }
  d_os << "END\n";
  d_os << "# end of topology file" << endl;
}

//OutTopology &OutTopology::operator<<(const gcore::Simulation &sim){
void OutTopology::write96(const gcore::System &sys, const gcore::GromosForceField &gff){
  
  // Title block
  d_os << "TITLE\n" << d_title << "\nEND\n";

  // TOPPHYSCON block
  d_os.precision(10);
  
  d_os << "TOPPHYSCON\n"
       << "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)\n"
       << gff.fpepsi() 
       << "\n# HBAR: Planck's constant HBAR = H/(2* PI)\n"
       << gff.hbar()
       <<"\nEND\n";
  
  // TOPVERSION block
  d_os << "TOPVERSION\n1.7\nEND\n";

  // ATOMTYPENAME block
  d_os << "ATOMTYPENAME\n"
       << "# NRATT: number of van der Waals atom types\n";
  int num=gff.numAtomTypeNames();
  d_os << num << "\n";
  d_os << "# TYPE: atom type names\n";
  for(int i=0;i<num;++i){
    if(i>0 &&!(i%10))d_os << "# " << i << "\n";
    d_os << gff.atomTypeName(i) << "\n";
  }
  d_os << "END\n";
  
  // RESNAME block
  d_os << "RESNAME\n"
       << "# NRAA2: number of residues in a solute molecule\n";
  num=0;
  for(int i=0;i<sys.numMolecules();++i)
    num+=sys.mol(i).topology().numRes();

  d_os << num << "\n"
       << "# AANM: residue names\n";
  for(int i=0, count=0;i<sys.numMolecules();++i)
    for(int j=0;j<sys.mol(i).topology().numRes();++j,++count){
      if(count>0 &&!(count%10))d_os << "# " << count << "\n";
      d_os << sys.mol(i).topology().resName(j) << "\n";
  }
  d_os << "END\n";
  
  // SOLUTEATOM block
  d_os << "SOLUTEATOM\n"
       << "#   NRP: number of solute atoms\n";
  num=0;
  for(int i=0;i<sys.numMolecules();++i)
    num+=sys.mol(i).numAtoms();

  d_os << setw(5) << num << "\n";
  d_os << "#  ATNM: atom number\n"
       << "#  MRES: residue number\n"
       << "#  PANM: atom name of solute atom\n"
       << "#   IAC: integer (van der Waals) atom type code\n"
       << "#  MASS: mass of solute atom\n"
       << "#    CG: charge of solute atom\n"
       << "#   CGC: charge group code (0 or 1)\n"
       << "#   INE: number of excluded atoms\n"
       << "# INE14: number of 1-4 interactions\n"
       << "# ATNM MRES PANM IAC     MASS       CG  CGC INE\n"
       << "#                                           INE14\n";

  for(int i=0, offatom=1, offres=1;i<sys.numMolecules();++i){
    for(int j=0;j<sys.mol(i).numAtoms();++j){
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(6) << offatom+j << ' '
	   << setw(4) << sys.mol(i).topology().resNum(j)+offres << ' '
	   << setw(4) << sys.mol(i).topology().atom(j).name()
	   << setw(4) << sys.mol(i).topology().atom(j).iac()+1
	   << setw(9) << sys.mol(i).topology().atom(j).mass()
	   << setw(9) << sys.mol(i).topology().atom(j).charge() 
	   << setw(3) << sys.mol(i).topology().atom(j).chargeGroup();
      // Exclusions
      d_os << setw(3) << sys.mol(i).topology().atom(j).exclusion().size();
      for(int k=0;k<sys.mol(i).topology().atom(j).exclusion().size();++k){
	if(k%6==0 && k!=0)
	  d_os << "\n"
	       << "                                            ";
	d_os << setw(6) << sys.mol(i).topology().atom(j).exclusion().atom(k)+offatom;
      }
      d_os << "\n"
	   << "                                           " 
	   << sys.mol(i).topology().atom(j).exclusion14().size();
      for(int k=0;k<sys.mol(i).topology().atom(j).exclusion14().size();++k){
	if(k%6==0 && k!=0)
	  d_os << "\n"
	       << "                                            ";
	d_os << setw(6) << sys.mol(i).topology().atom(j).exclusion14().atom(k)+offatom;
      }
      
      d_os << "\n";
    }
    offres+=sys.mol(i).topology().numRes();
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // BONDTYPE block

  d_os << "BONDTYPE\n"
       << "#  NBTY: number of covalent bond types\n";
  num=gff.numBondTypes();
  
  d_os << num << "\n"
       << "#  CB: force constant\n"
       << "#  B0: bond length at minimum energy\n"
       << "#         CB          B0\n";

  for (int i=0;i<num;++i){
    if(i>0 &&!(i%10)) d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.bondType(i).fc()
	 << setw(12) << gff.bondType(i).b0() << "\n";
  }
  d_os << "END\n";
  
  // BONDH block
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH())
	++num;
    }
  }
  d_os << "BONDH\n"
       << "#  NBONH: number of bonds involving H atoms in solute\n"
       << num << "\n"
       << "#  IBH, JBH: atom sequence numbers of atoms forming a bond\n"
       << "#  ICBH: bond type code\n"
       << "#   IBH    JBH ICBH\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";
  
  d_os << "BOND\n"
       << "#  NBON: number of bonds NOT involving H atoms in solute\n";
  num=0;
  for(int i=0; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH()) 
	++num;
    }
  }
  d_os << num << "\n"
       << "#  IB, JB: atom sequence numbers of atoms forming a bond\n"
       << "#  ICB: bond type code\n"
       << "#    IB     JB  ICB\n";
  
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    BondIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH()) {
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      } 
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // BONDANGLETYPE block
  num=gff.numAngleTypes();
  d_os << "BONDANGLETYPE\n"
       << "#  NTTY: number of bond angle types\n"
       << num << "\n"
       << "#  CT: force constant\n"
       << "#  T0: bond angle at minimum energy in degrees\n"
       << "#         CT          T0\n";

  for (int i=0;i<num;++i){
    if(i>0 &&!(i%10))d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.angleType(i).fc()
	 << setw(12) << gff.angleType(i).t0() << "\n";
  }
  d_os << "END\n";

  // BONDANGLEH & BONDANGLE block
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
  d_os << "BONDANGLEH\n"
       << "#  NTHEH: number of bond angles involving H atoms in solute\n"
       << num << "\n"
       << "#  ITH, JTH, KTH: atom sequence numbers\n"
       << "#    of atoms forming a bond angle in solute\n"
       << "#  ICTH: bond angle type code\n"
       << "#   ITH    JTH    KTH ICTH\n";
  
  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    AngleIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH()){
	
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2]+ offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

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
  d_os << "BONDANGLE\n"
       << "#  NTHE: number of bond angles NOT\n"
       << "#        involving H atoms in solute\n"
       << num << "\n"
       << "#  IT, JT, KT: atom sequence numbers of atoms\n"
       << "#     forming a bond angle\n"
       << "#  ICT: bond angle type code\n"
       << "#    IT     JT     KT  ICT\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    AngleIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2] +offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";
      
  // IMPDIHEDRALTYPE block
  num=gff.numImproperTypes();
  d_os << "IMPDIHEDRALTYPE\n"
       << "#  NQTY: number of improper dihedrals\n"
       << num << "\n"
       << "#  CQ: force constant of improper dihedral per degrees square\n"
       << "#  Q0: improper dihedral angle at minimum energy in degrees\n"
       << "#         CQ          Q0\n";

  for (int i=0;i<num;++i){
    if(i>0 &&!(i%10))d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.improperType(i).fc()
	 << setw(12) << gff.improperType(i).q0() << "\n";
  }
  d_os << "END\n";

  // IMPDIHEDRALH & IMPDIHEDRAL block
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
  d_os << "IMPDIHEDRALH\n"
       << "#  NQHIH: number of improper dihedrals\n"
       << "#         involving H atoms in the solute\n"
       << num << "\n"
       << "#  IQH,JQH,KQH,LQH: atom sequence numbers\n"
       << "#     of atoms forming an improper dihedral\n"
       << "#  ICQH: improper dihedral type code\n"
       << "#   IQH    JQH    KQH    LQH ICQH\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    ImproperIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH() ||
	 sys.mol(i).topology().atom(bit()[3]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2]+ offatom
	     << setw(7) << bit()[3]+ offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

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
  d_os << "IMPDIHEDRAL\n"
       << "#  NQHI: number of improper dihedrals NOT\n"
       << "#    involving H atoms in solute\n"
       << num << "\n"
       << "#  IQ,JQ,KQ,LQ: atom sequence numbers of atoms\n"
       << "#    forming an improper dihedral\n"
       << "#  ICQ: improper dihedral type code\n"
       << "#    IQ     JQ     KQ     LQ  ICQ\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    ImproperIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH() &&
	 !sys.mol(i).topology().atom(bit()[3]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2] +offatom
	     << setw(7) << bit()[3] + offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // DIHEDRALTYPE block
  num=gff.numDihedralTypes();

  d_os << "DIHEDRALTYPE\n"
       << "#  NPTY: number of dihedral types\n"
       << num << "\n"
       << "#  CP: force constant\n"
       << "#  PD: cosine of the phase shift\n"
       << "#  NP: multiplicity\n"
       << "#       CP        PD  NP\n";

  for (int i=0;i<num;++i){
    if(i>0 &&!(i%10))d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os << setw(10) << gff.dihedralType(i).fc()
	 << setw(10) << gff.dihedralType(i).pd() 
	 << setw(4)<< gff.dihedralType(i).np() << "\n";
  }
  d_os << "END\n";

  // DIHEDRALH & DIHEDRAL block
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
  d_os << "DIHEDRALH\n"
       << "#  NPHIH: number of dihedrals involving H atoms in solute\n"
       << num << "\n"
       << "#  IPH, JPH, KPH, LPH: atom sequence numbers\n"
       << "#    of atoms forming a dihedral\n"
       << "#  ICPH: dihedral type code\n"
       << "#   IPH    JPH    KPH    LPH ICPH\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    DihedralIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(sys.mol(i).topology().atom(bit()[0]).isH() ||
	 sys.mol(i).topology().atom(bit()[1]).isH() ||
	 sys.mol(i).topology().atom(bit()[2]).isH() ||
	 sys.mol(i).topology().atom(bit()[3]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2]+ offatom
	     << setw(7) << bit()[3]+ offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";

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
  d_os << "DIHEDRAL\n"
       << "#  NPHI: number of dihedrals NOT involving H atoms in solute\n"
       << num << "\n"
       << "#  IP, JP, KP, LP: atom sequence numbers\n"
       << "#     of atoms forming a dihedral\n"
       << "#  ICP: dihedral type code\n"
       << "#    IP     JP     KP     LP  ICP\n";

  for(int i=0, offatom=1; i<sys.numMolecules(); ++i){
    DihedralIterator bit(sys.mol(i).topology());
    for(int count=0;bit;++bit){
      if(!sys.mol(i).topology().atom(bit()[0]).isH() &&
	 !sys.mol(i).topology().atom(bit()[1]).isH() &&
	 !sys.mol(i).topology().atom(bit()[2]).isH() &&
	 !sys.mol(i).topology().atom(bit()[3]).isH()){
	if(count>0 &&!(count%10))d_os << "# " << count << "\n";
	d_os << setw(7) << bit()[0] +offatom
	     << setw(7) << bit()[1]+offatom
	     << setw(7) << bit()[2] +offatom
	     << setw(7) << bit()[3] + offatom
	     << setw(5) << bit().type()+1 << "\n";
	++count;
      }
    }
    offatom+=sys.mol(i).numAtoms();
  }
  d_os << "END\n";
  
  // LJPARAMETERS block
  num=gff.numLJTypes();
  d_os << "LJPARAMETERS\n"
       << "#  NRATT2: number of LJ interaction types = NRATT*(NRATT+1)/2\n"
       << num << "\n"
       << "#  IAC,JAC: integer (van der Waals) atom type code\n"
       << "#  C12: r**(-12) term in nonbonded interactions\n"
       << "#   C6: r**(-6) term in nonbonded interactions\n"
       << "# CS12: r**(-12) term in 1-4 nonbonded interactions\n"
       << "#  CS6: r**(-6) term in 1-4 nonbonded interactions\n"
       << "# IAC  JAC           C12            C6          CS12           CS6\n";

  for (int i=0;i<gff.numAtomTypeNames();++i){
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
	   << setw(14) << lj.cs6() << "\n";
    }
    d_os << "#\n";
  }
  d_os << "END\n";

  //SOLVENTATOM block
  d_os << "SOLVENTATOM\n"
       << "#  NRAM: number of atoms per solvent molecule\n"
       << sys.sol(0).topology().numAtoms() << "\n"
       << "#     I: solvent atom sequence number\n"
       << "#  IACS: integer (van der Waals) atom type code\n"
       << "#  ANMS: atom name of solvent atom\n"
       << "#  MASS: mass of solvent atom\n"
       << "#   CGS: charge of solvent atom\n"
       << "#  I  ANMS IACS      MASS        CGS\n";

  for(int j=0;j<sys.sol(0).topology().numAtoms();++j){
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os << setw(4) << j+1 << " "
	 << setw(5) << sys.sol(0).topology().atom(j).name()
	 << setw(4) << sys.sol(0).topology().atom(j).iac()+1
	 << setw(11) << sys.sol(0).topology().atom(j).mass()
	 << setw(11) << sys.sol(0).topology().atom(j).charge();
    d_os << "\n";
  }
  d_os << "END\n";

  //SOLVENTCONSTR bock
  num=0;
  ConstraintIterator dit(sys.sol(0).topology());
  for(;dit;++dit) ++num;
  d_os << "SOLVENTCONSTR\n"
       << "#  NCONS: number of constraints\n"
       << num << "\n"
       << "#  ICONS, JCONS: atom sequence numbers forming constraint\n"
       << "#   CONS constraint length\n"
       << "#ICONS JCONS         CONS\n";

  ConstraintIterator cit(sys.sol(0).topology());

  for(;cit;++cit){
    d_os.precision(7);
    d_os.setf(ios::fixed, ios::floatfield);

    d_os << setw(5) << cit()[0] + 1
         << setw(5) << cit()[1] + 1
         << setw(15) << cit().dist() << "\n";
  }
  d_os << "END\n";
  d_os << "# end of topology file" << endl;
}

