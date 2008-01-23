/**
 * @file make_pt_top.cc
 * Creates perturbation topologies from molcular topologies
 */

/**
 * @page programs Program Documentation
 *
 * @anchor make_pt_top
 * @section make_pt_top Create a perturbation topology from two molecular topologies
 * @author @ref ns
 * @date 22.01.08
 * Program make_pt_top takes two molecular topologies and writes the 
 * differences in the perturbation topology format. Both topologies must contain
 * the same number of solute atoms. The softness parameters @f$\alpha_{LJ}@f$ 
 * @f$\alpha_{CRF}@f$ can be specified. 
 * 
 * Atoms, atom pairs, bonds, bond angles, 
 * improper dihedrals and dihedrals contained in the \@select atom specifier are
 * kept even if the parameters in both topologies are the same. Atoms, atom 
 * pairs, bonds, bond angles and dihedrals contained in the \@reject atom
 * specifier are deleted from the output.
 * The resulting perturbation topology file (in PERTOPO03 format) is written out
 * to the standard output.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topoA</td><td>&lt;molecular topology file for state A&gt; </td></tr>
 * <tr><td> \@topoB</td><td>&lt;molecular topology file for state B&gt; </td></tr>
 * <tr><td> \@softpar</td><td>&lt;softness parameters (@f$\alpha_{LJ}@f$ @f$\alpha_{CRF}@f$)&gt; </td></tr>
 * <tr><td> \@select</td><td>&lt;@ref AtomSpecifier "AtomSpecifier": atoms to include&gt; </td></tr>
 * <tr><td> \@reject</td><td>&lt;@ref AtomSpecifier "AtomSpecifier": atoms to reject&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  make_top
    @topoA    exA.top
    @topoB    exB.top
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/utils/AtomSpecifier.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

#include "../src/utils/make_pt_top.h"

int main(int argc, char *argv[]){

  char *knowns[] = {"topoA", "topoB", "softpar", "select", "reject"};
  int nknowns = 5;
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topoA      <topology for state A>\n";
  usage += "\t@topoB      <topology for state B>\n";
  usage += "\t[@softpar   <alpha_lj alpha_crf (default 1.51 0.5)>]\n"; 
  usage += "\t[@select    <Atomspecifier atoms to keep>]\n";
  usage += "\t[@reject    <Atomspecifier atoms not to keep>]\n";
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    InTopology itA(args["topoA"]), itB(args["topoB"]);
    System sysA(itA.system()), sysB(itB.system());
    LinearTopology topA(sysA), topB(sysB);
    
    // abort if the topologies do not contain the same number of atoms.
    // mapping is risky and too difficult.
    if (topA.atoms().size() != topB.atoms().size())
      throw gromos::Exception("make_pt_top", "Both topologies need to have "
              "the same number of atoms.");
    
    double alpha_lj = 1.51, alpha_crf = 0.5;
    // get the softness parameters
    if (args.count("softpar") > 0) {
      if (args.count("softpar") != 2) {
        throw gromos::Exception("make_pt_top", "Both softness parameters "
                "(alpha_LJ and alpha_CRF) have to be given.");
      }
      Arguments::const_iterator iter=args.lower_bound("softpar");
      istringstream arg1(iter->second), arg2((++iter)->second);
      if (!(arg1 >> alpha_lj) || !(arg2 >> alpha_crf)) {
        throw gromos::Exception("make_pt_top", "Make sure the softness "
                "parameters are floating point numbers.");
      }
    }
    
    // read in the atom list the has to be kept definitely
    utils::AtomSpecifier select(sysA);
    {
      Arguments::const_iterator iter=args.lower_bound("select"),
	to=args.upper_bound("select");
      for(;iter!=to; ++iter)
	select.addSpecifier(iter->second);
    }    
    // read in the atom list for atoms to throw away for sure
    utils::AtomSpecifier reject(sysA);
    {
      Arguments::const_iterator iter=args.lower_bound("reject"),
	to=args.upper_bound("reject");
      for(;iter!=to; ++iter)
	reject.addSpecifier(iter->second);
    }
    // check that there are no doubles in ls and rej
    for(int i=0; i<reject.size(); i++)
      if( select.findAtom(reject.mol(i), reject.atom(i))!=-1 )
	throw gromos::Exception("make_pt_top", "select and reject show overlap");
    
    cout << "TITLE" << endl << argv[0]
         << " generated perturbation topology." << endl
         << "State A: " << args["topoA"] << endl
         << "State B: " << args["topoB"] << endl
         << "END" << endl
         << "SCALEDINTERACTIONS" << endl << "0" << endl << "END" << endl;
    
    vector<PertAtom> pas;
    set<PertAtomPair> paps; // set, because they may be inserted twice.
    
    // loop over all atoms
    for(unsigned int i = 0; i < topA.atoms().size(); i++) {
      AtomTopology& atomA = topA.atoms()[i], atomB = topB.atoms()[i];
      
      if (((atomA.iac() != atomB.iac() ||
            atomA.mass() != atomB.mass() ||
            atomA.charge() != atomB.charge()) &&
           !containsGromosNum(reject, i)) || containsGromosNum(select, i)) {
        if (atomA.name() != atomB.name()) {
          // this is not bad but worth a warning
          cerr << "Warning: atom " << i+1 << " name mismatch (" <<
                  atomA.name() << " != " << atomB.name() << ")" << endl;        
        }
        
        // store data in a PertAtom struct
        PertAtom pa(alpha_lj, alpha_crf);
        pa.gromosNum = i;
        pa.name = atomA.name();
        pa.resNum = topA.resMap()[i];
        pa.iac[0] = atomA.iac();
        pa.iac[1] = atomB.iac();
        pa.mass[0] = atomA.mass();
        pa.mass[1] = atomB.mass();
        pa.charge[0] = atomA.charge();
        pa.charge[1] = atomB.charge();
        
        pas.push_back(pa);
      }
      
      // exclusions are difficult:
      // check wheter all atoms excluded by a certain atom in topology A
      // are also excluded in topology B and vice versa. 
      // don't forget 1-4 exclusions.
      const Exclusion& exA = atomA.exclusion(), exB = atomB.exclusion(),
                 ex14A = atomA.exclusion14(), ex14B = atomB.exclusion14();
      // exclusions of atomA also in atomB?
      for(int j = 0; j < exA.size(); j++) {
        unsigned stateA = 0; // excluded
        unsigned stateB = 1; // normal interaction
        if (exclusionContains(exB, exA.atom(j)))
          stateB = 0;
        else if (exclusionContains(ex14B, exA.atom(j)))
          stateB = 2;
        
        if (stateA != stateB || containsGromosNum(select, i) || containsGromosNum(select, exA.atom(j))) {
          PertAtomPair pap;
          pap.gromosNum[0] = i;
          pap.gromosNum[1] = exA.atom(j);
          pap.interactionType[0] = stateA;
          pap.interactionType[1] = stateB;
          if (!containsGromosNum(reject, pap.gromosNum[0]) &&
              !containsGromosNum(reject, pap.gromosNum[1]))
            paps.insert(pap);
        }                 
      }
      // exclusions14 of atomA also in atomB?
      for(int j = 0; j < ex14A.size(); j++) {
        unsigned stateA = 2; // 1-4 excluded
        unsigned stateB = 1; // normal interaction
        if (exclusionContains(exB, ex14A.atom(j)))
          stateB = 0;
        else if (exclusionContains(ex14B, ex14A.atom(j)))
          stateB = 2;
        
        if (stateA != stateB || containsGromosNum(select, i) || containsGromosNum(select, ex14A.atom(j))) {
          PertAtomPair pap;
          pap.gromosNum[0] = i;
          pap.gromosNum[1] = ex14A.atom(j);
          pap.interactionType[0] = stateA;
          pap.interactionType[1] = stateB;
          if (!containsGromosNum(reject, pap.gromosNum[0]) &&
              !containsGromosNum(reject, pap.gromosNum[1]))
            paps.insert(pap);
        }                 
      }
      // and the same the other way round?
      // exclusions of atomB also in atomA?
      for(int j = 0; j < exB.size(); j++) {
        unsigned stateB = 0; // excluded
        unsigned stateA = 1; // normal interaction
        if (exclusionContains(exA, exB.atom(j)))
          stateA = 0;
        else if (exclusionContains(ex14A, exB.atom(j)))
          stateA = 2;
        
        if (stateA != stateB || containsGromosNum(select, i) || containsGromosNum(select, exB.atom(j))) {
          PertAtomPair pap;
          pap.gromosNum[0] = i;
          pap.gromosNum[1] = exB.atom(j);
          pap.interactionType[0] = stateA;
          pap.interactionType[1] = stateB;
          if (!containsGromosNum(reject, pap.gromosNum[0]) &&
              !containsGromosNum(reject, pap.gromosNum[1]))
            paps.insert(pap);
        }                 
      }
      // exclusions14 of atomB also in atomA?
      for(int j = 0; j < ex14B.size(); j++) {
        unsigned stateB = 2; // 1-4 excluded
        unsigned stateA = 1; // normal interaction
        if (exclusionContains(exA, ex14B.atom(j)))
          stateA = 0;
        else if (exclusionContains(ex14A, ex14B.atom(j)))
          stateA = 2;
        
        if (stateA != stateB || containsGromosNum(select, i) || containsGromosNum(select, ex14B.atom(j))) {
          PertAtomPair pap;
          pap.gromosNum[0] = i;
          pap.gromosNum[1] = ex14B.atom(j);
          pap.interactionType[0] = stateA;
          pap.interactionType[1] = stateB;
          if (!containsGromosNum(reject, pap.gromosNum[0]) &&
              !containsGromosNum(reject, pap.gromosNum[1]))
            paps.insert(pap);
        }                 
      }
    }
    
    // write out PERTATOM03 block
    cout << "PERTATOM03" << endl
         << "# number of perturbed atoms" << endl
         << pas.size() << endl
         << "#   NR RES NAME IAC(A) MASS(A) CHARGE(A) IAC(B) MASS(B) "
            "CHARGE(B) ALJ   ACRF" << endl;
    // here we can use the copy algorithm with ostream_iterator because
    // all Pert* structs have an overloaded << operator.
    copy(pas.begin(), pas.end(), ostream_iterator<PertAtom>(cout));
    cout << "END" << endl;
    
    // write out PERTATOMPAIR03
    cout << "PERTATOMPAIR03" << endl
         << "# number of perturbed atom pairs" << endl
         << paps.size() << endl
         << "#    i     j i(A) i(B)" << endl;
    copy(paps.begin(), paps.end(), ostream_iterator<PertAtomPair>(cout));
    cout << "END" << endl;
         
    vector<PertBond> pbs;
    
    if (topA.bonds().size() != topB.bonds().size()) {
      cerr << "Warning: the topologies don't have the same "
              "number of bonds!" << endl
           << "         Taking bonds of state A and searching "
              "for changes in state B." << endl;
    }
    
    // loop over all bonds
    set<Bond>::const_iterator itBondA = topA.bonds().begin(), 
            toBondA = topA.bonds().end();
    
    for(; itBondA != toBondA; ++itBondA) {
      const Bond& bondA = *itBondA;
      bool found = false;
      
      // search for bondA in topology B -> loop over bonds in topology B
      set<Bond>::const_iterator itBondB = topB.bonds().begin(), 
            toBondB = topB.bonds().end();
      
      for(; itBondB != toBondB; ++itBondB) {
        const Bond& bondB = *itBondB;
        if ((bondA[0] == bondB[0] && bondA[1] == bondB[1]) ||
            (bondA[0] == bondB[1] && bondA[1] == bondB[0])) {
          found = true;
          if ((bondA.type() != bondB.type() && !containsGromosNum(reject, bondA)) ||
               containsGromosNum(select, bondA)) {
            // store the data
            PertBond pb;
            pb.gromosNum[0] = bondA[0];
            pb.gromosNum[1] = bondA[1];
            pb.type[0] = bondA.type();
            pb.type[1] = bondB.type();
            pbs.push_back(pb);
          }
          break;
        }
      }
      if (!found) {
        cerr << "Warning: could not find bondA " << bondA[0]+1 << "-"
             << bondA[1]+1 << endl;
      }
    }
    
    // write out PERTBOND03 block
    cout << "PERTBOND03" << endl
         << "# number of perturbed bonds" << endl
         << pbs.size() << endl
         << "#    i     j t(A) t(B)" << endl;
    copy(pbs.begin(), pbs.end(), ostream_iterator<PertBond>(cout));
    cout << "END" << endl;
    
    vector<PertAngle> pangs;
    
    if (topA.angles().size() != topB.angles().size()) {
      cerr << "Warning: the topologies don't have the same "
              "number of angles!" << endl
           << "         Taking bond angles of state A and searching "
              "for changes in state B." << endl;
    }
    
    // loop over all angles
    set<Angle>::const_iterator itAngleA = topA.angles().begin(), 
            toAngleA = topA.angles().end();
    
    for(; itAngleA != toAngleA; ++itAngleA) {
      const Angle& angleA = *itAngleA;
      bool found = false;
      
      // search for angleA in topology B
      set<Angle>::const_iterator itAngleB = topB.angles().begin(), 
            toAngleB = topB.angles().end();
      
      for(; itAngleB != toAngleB; ++itAngleB) {
        const Angle& angleB = *itAngleB;
        if (angleA[0] == angleB[0] &&
            angleA[1] == angleB[1] &&
            angleA[2] == angleB[2]) {
          found = true;
          if ((angleA.type() != angleB.type() && !containsGromosNum(reject, angleA)) ||
              containsGromosNum(select, angleA)) {
            // store the data
            PertAngle pa;
            pa.gromosNum[0] = angleA[0];
            pa.gromosNum[1] = angleA[1];
            pa.gromosNum[2] = angleA[2];
            pa.type[0] = angleA.type();
            pa.type[1] = angleB.type();
            pangs.push_back(pa);
          }
          break;
        }
      }
      if (!found) {
        cerr << "Warning: could not find angleA " << angleA[0]+1 << "-"
             << angleA[1]+1 << "-" << angleA[2]+1 << endl;
      }
    }
    
    // write out PERANGLE03
    cout << "PERTANGLE03" << endl
         << "# number of perturbed bond angles" << endl
         << pangs.size() << endl
         << "#    i     j     k t(A) t(B)" << endl;
    copy(pangs.begin(), pangs.end(), ostream_iterator<PertAngle>(cout));
    cout << "END" << endl;
    
    vector<PertImproper> pimps;
  
    if (topA.impropers().size() != topB.impropers().size()) {
      cerr << "Warning: the topologies don't have the same "
              "number of improper dihedrals!" << endl
           << "         Taking improper dihedrals of state A and searching "
              "for changes in state B." << endl;
    }
    
    // loop over impropers
    set<Improper>::const_iterator itImproperA = topA.impropers().begin(), 
            toImproperA = topA.impropers().end();
    
    for(; itImproperA != toImproperA; ++itImproperA) {
      const Improper& improperA = *itImproperA;
      bool found = false;
      
      // search for improperA in topology B
      set<Improper>::const_iterator itImproperB = topB.impropers().begin(), 
            toImproperB = topB.impropers().end();
      
      for(; itImproperB != toImproperB; ++itImproperB) {
        const Improper& improperB = *itImproperB;
        if (improperA[0] == improperB[0] &&
            improperA[1] == improperB[1] &&
            improperA[2] == improperB[2] &&
            improperA[3] == improperB[3]) {
          found = true;
          if ((improperA.type() != improperB.type() && !containsGromosNum(reject, improperA)) ||
              containsGromosNum(select, improperA)) {
            // store the data
            PertImproper pi;
            pi.gromosNum[0] = improperA[0];
            pi.gromosNum[1] = improperA[1];
            pi.gromosNum[2] = improperA[2];
            pi.gromosNum[3] = improperA[3];
            pi.type[0] = improperA.type();
            pi.type[1] = improperB.type();
            pimps.push_back(pi);
          }
          break;
        }
      }
      if (!found) {
        cerr << "Warning: could not find improperA " << improperA[0]+1 << "-"
             << improperA[1]+1 << "-" << improperA[2]+1 << endl;
      }
    }
    
    // write PERTIMPDIHEDRAL03 block
    cout << "PERTIMPDIHEDRAL03" << endl
         << "# number of perturbed improper dihedrals" << endl
         << pimps.size() << endl
         << "#    i     j     k     l t(A) t(B)" << endl;
    copy(pimps.begin(), pimps.end(), ostream_iterator<PertImproper>(cout));
    cout << "END" << endl;
    
    vector<PertDihedral> pds;
  
    if (topA.dihedrals().size() != topB.dihedrals().size()) {
      cerr << "Warning: the topologies don't have the same "
              "number of dihedral dihedrals!" << endl
           << "         Taking dihedrals of state A and searching "
              "for changes in state B." << endl;
    }
    
    // loop over dihedrals
    set<Dihedral>::const_iterator itDihedralA = topA.dihedrals().begin(), 
            toDihedralA = topA.dihedrals().end();
    
    for(; itDihedralA != toDihedralA; ++itDihedralA) {
      const Dihedral& dihedralA = *itDihedralA;
      bool found = false;
      
      // search for dihedralA in topology B
      
      set<Dihedral>::const_iterator itDihedralB = topB.dihedrals().begin(), 
            toDihedralB = topB.dihedrals().end();
      
      for(; itDihedralB != toDihedralB; ++itDihedralB) {
        const Dihedral& dihedralB = *itDihedralB;
        if (dihedralA[0] == dihedralB[0] &&
            dihedralA[1] == dihedralB[1] &&
            dihedralA[2] == dihedralB[2] &&
            dihedralA[3] == dihedralB[3]) {
          found = true;
          if ((dihedralA.type() != dihedralB.type() && !containsGromosNum(reject, dihedralA)) ||
              containsGromosNum(select, dihedralA)) {
            // store data
            PertDihedral pd;
            pd.gromosNum[0] = dihedralA[0];
            pd.gromosNum[1] = dihedralA[1];
            pd.gromosNum[2] = dihedralA[2];
            pd.gromosNum[3] = dihedralA[3];
            pd.type[0] = dihedralA.type();
            pd.type[1] = dihedralB.type();
            pds.push_back(pd);
          }
          break;
        }
      }
      if (!found) {
        cerr << "Warning: could not find dihedralA " << dihedralA[0]+1 << "-"
             << dihedralA[1]+1 << "-" << dihedralA[2]+1 << endl;
      }
    }
    // wirte PERTDIHEDRAL03 block
    cout << "PERTDIHEDRAL03" << endl
         << "# number of perturbed dihedrals" << endl
         << pds.size() << endl
         << "#    i     j     k     l t(A) t(B)" << endl;
    copy(pds.begin(), pds.end(), ostream_iterator<PertDihedral>(cout));
    cout << "END" << endl;
    
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}
