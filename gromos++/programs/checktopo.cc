//checktopo reads in a topology and performs some simple checks on it
//          if given a building block file, it can check your topology for
//          consistency. If a coordinate file is given, it will also write out
//          the energy for all bonded interaction. (This feature will go to 
//          shake_analysis in the future)

#include <cassert>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gio/InParameter.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Energy.h"
#include "../src/utils/CheckTopo.h"
#include "../src/utils/FfExpert.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

  
int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "coord", "build", "param"};
  int nknowns =5;

  string usage = argv[0];
  usage += "\n\t@topo     <topology>\n";
  usage += "\t[@coord   <coordinate file>\n";
  usage += "\t[@pbc     <boundary type> <gather method>\n";
  usage += "\t[@build   <building block file for consistency check>\n";
  usage += "\t[@param   <parameter file for consistency check>\n";
  
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  //  read topology
  InTopology it(args["topo"]);
  System sys(it.system());
  int nummol=sys.numMolecules();
  GromosForceField gff(it.forceField());
  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);
  Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
  
  // declare the energy class
  Energy en(sys, gff, *pbc);
  
  // set properties
  PropertyContainer props(sys);
  int numbonds[nummol];
  int numangles[nummol];
  int numimp[nummol];
  int numdih[nummol];
  int maxatomtype=0, maxbondtype=0, maxangletype=0, maximptype=0, maxdihtype=0;
  
  double totcharge[nummol];

  // loop over all bonds
  for(int m=0; m<nummol; m++){
    numbonds[m]=0;
    BondIterator bi(sys.mol(m).topology());
    for(;bi;++bi){
      ostringstream os;
      os << "d%" << m+1 << ":" << bi()[0]+1 << "," << bi()[1]+1;
      props.addSpecifier(os.str());
      numbonds[m]++;
      if(bi().type()>maxbondtype) maxbondtype=bi().type();
    }
  }
  // loop over all angles
  for(int m=0; m<nummol; m++){
    numangles[m]=0;
    AngleIterator ai(sys.mol(m).topology());
    for(;ai;++ai){
      ostringstream os;
      os << "a%" << m+1 << ":" << ai()[0]+1 << "," << ai()[1]+1 << "," 
	 << ai()[2]+1;
      props.addSpecifier(os.str());
      numangles[m]++;
      if(ai().type()>maxangletype) maxangletype=ai().type();
    }
  }
  // loop over all impropers
  for(int m=0; m<nummol; m++){
    numimp[m]=0;
    ImproperIterator ii(sys.mol(m).topology());
    for(;ii;++ii){
      ostringstream os;
      os << "t%" << m+1 << ":" << ii()[0]+1 << "," << ii()[1]+1 << "," 
	 << ii()[2]+1 << "," << ii()[3]+1;
      props.addSpecifier(os.str());
      numimp[m]++;
      if(ii().type()>maximptype) maximptype=ii().type();
    }
  }
  // loop over all dihedrals
  for(int m=0; m<nummol; m++){
    numdih[m]=0;
    DihedralIterator di(sys.mol(m).topology());
    for(;di;++di){
      ostringstream os;
      os << "t%" << m+1 << ":" << di()[0]+1 << "," << di()[1]+1 << "," 
	 << di()[2]+1 << "," << di()[3]+1;
      props.addSpecifier(os.str());
      numdih[m]++;
      if(di().type()>maxdihtype) maxdihtype=di().type();
    }
  }

  // calculate the total charge
  for(int m=0; m<nummol; m++){
    totcharge[m]=0.0;
    for(int a=0; a<sys.mol(m).numAtoms(); a++){
      totcharge[m]+=sys.mol(m).topology().atom(a).charge();
      if(sys.mol(m).topology().atom(a).iac() > maxatomtype)
	maxatomtype=sys.mol(m).topology().atom(a).iac();
    }
  }
  
  if(args.count("coord")>0){
    
    // now, we are done preparing everything the real program starts here
    // calculate the values of all the properties
    props.calc();
    
    // parse them into the energy class
    en.setProperties(props);
    
    // define input coordinate
    InG96 ic;
    
    // open file
    ic.open(args["coord"]);
    ic.select("SOLUTE");
    
    // read coordinates and gather (gather method does not really matter)
    ic >> sys;
    (*pbc.*gathmethod)();
    
    // calculate the energies
    en.calc();
  }
  
  // That was it, now for the output
  int tnumbonds=0, tnumangles=0, tnumimp=0, tnumdih=0;
  double ttotcharge=0.0;
  double chrg_precision=0.0001;
  cout << "Topology contains " << sys.numMolecules() << " molecule";
  if(sys.numMolecules()>1) cout << "s";
  cout << ":" << endl << endl;
  cout << setw(10) << "molecule"
       << setw(12) << "# atoms"
       << setw(12) << "# bonds"
       << setw(12) << "# angles"
       << setw(12) << "# impropers"
       << setw(12) << "# dihedrals" 
       << setw(12) << "tot charge" << endl;
  for(int m=0; m< nummol; m++){
    if(fabs(totcharge[m])<chrg_precision) totcharge[m]=0.0;
    cout  << setw(10) << m+1
	  << setw(12) << sys.mol(m).topology().numAtoms()
	  << setw(12) << numbonds[m]
	  << setw(12) << numangles[m]
	  << setw(12) << numimp[m]
	  << setw(12) << numdih[m] 
	  << setw(12) << totcharge[m] << endl;
    tnumbonds+=numbonds[m];
    tnumangles+=numangles[m];
    tnumimp+=numimp[m];
    tnumdih+=numdih[m];
    ttotcharge+=totcharge[m];
    
  }
  
  // do some tests on the topology
  cout << endl 
       << "Performing some basic checks on the charge groups, exclusions, "
       << "bonds, angles and improper dihedrals..." << endl;
  int error = 0;
  
  // loop over the molecules
  for(int m=0; m< sys.numMolecules(); m++){

    // define a topology check
    utils::CheckTopo ct(sys,m);
    ct.checkAll();
    error += ct.numErrors();
    if(ct.numErrors()){
      cout << "--------------------" << endl;
      cout << "In Molecule " << m+1 << ":" << endl << endl;
      for(int i=0; i< ct.numErrors(); i++)
	cout << ct.error(i) << endl << endl;
      cout << "--------------------" << endl;
    }
  }
  // check whether all types are defined
  ostringstream par;
  int parerr=0;
  if(maxatomtype >= gff.numAtomTypeNames()){
    par << "Higher atom type used than defined: " << maxatomtype + 1
	<< " > " << gff.numAtomTypeNames() << endl << endl;
    ++parerr;
  }
  if(maxbondtype >= gff.numBondTypes()){
    par << "Higher bond type used than defined: " << maxbondtype + 1
	<< " > " << gff.numBondTypes() << endl << endl;
    ++parerr;
  }
  if(maxangletype >= gff.numAngleTypes()){
    par << "Higher angle type used than defined: " << maxangletype + 1
	<< " > " << gff.numAngleTypes() << endl << endl;
    ++parerr;
  }
  if(maximptype >= gff.numImproperTypes()){
    par << "Higher improper dihedral type used than defined: " << maximptype+1 
	<< " > " << gff.numImproperTypes() << endl << endl;
    ++parerr;
  }
  if(maxdihtype >= gff.numDihedralTypes()){
    par << "Higher dihedral type used than defined: " << maxdihtype + 1
	<< " > " << gff.numDihedralTypes() << endl << endl;
    ++parerr;
  }
  error += parerr;
  if(parerr){
    cout << "--------------------" << endl;
    cout << "In Parameters:" << endl << endl;
    cout << par.str();
    cout << "--------------------" << endl;
  }
  if(error==0)
    cout << "ok" << endl;

  // possibly perform consistency check with the building block
  if(args.count("build")>0){
    if(args.count("param")<=0)
      throw gromos::Exception("checktopo", 
			      "For consistency check, give both @buid and "
			      "@param input flags");
    cout << "\n\nComparing parameters with other building blocks for "
	 << "consistency\n"
	 << "Building block file: " << args["build"]
	 << "\nParameter file: " << args["param"] <<  endl;
    
    gio::InBuildingBlock ibb(args["build"]);
    gcore::BuildingBlock mtb(ibb.building());
    gio::InParameter ip(args["param"]);
    gcore::GromosForceField gffc(ip.forceField());
    
    utils::FfExpert ffexp(mtb);
    vector<utils::FfExpert::counter> v;
    
    // loop over all atoms
    for(int m=0; m < nummol; ++m){
      cout << "\n--------------------" << endl;
      cout << "In Molecule " << m+1 << ":\n\n";
      
      cout << endl << sys.mol(m).numAtoms() << " ATOMS :" << endl << endl;
      cout << setw(8) << "atom-"
	   << setw(8) << "atom-"
	   << setw(6) << "IAC"
	   << setw(9) << "mass"
	   << setw(9) << "charge"
	   << " charge group    consistency" << endl;
      cout << setw(8) << "number"
	   << setw(8) << "name" << endl;
      
      for(int a=0; a < sys.mol(m).numAtoms(); ++a){
	bool inconsistency=false;
	ostringstream wrn;
	
	// to prepare, write the atom information to the warning
	wrn << setw(8) << a+1 << ' '
	    << setw(7) << sys.mol(m).topology().atom(a).name()
	    << setw(6) << sys.mol(m).topology().atom(a).iac()+1
	    << setw(9) << sys.mol(m).topology().atom(a).mass()
	    << setw(9) << sys.mol(m).topology().atom(a).charge()
	    << setw(4) << sys.mol(m).topology().atom(a).chargeGroup();
	
	
	// check the name with the IAC
	ffexp.name2iac(sys.mol(m).topology().atom(a).name().substr(0,1), v);
	
	bool found=false;
	for(unsigned int i=0; i< v.size(); ++i){
	  if(v[i].type == sys.mol(m).topology().atom(a).iac()) found=true;
	}
	if(!found){
	  if(!inconsistency)
	    wrn << "             Inconsistency with building block file:\n";
	  wrn << "\t\tNo atoms found with name " 
	      << sys.mol(m).topology().atom(a).name() << " and IAC "
	      << sys.mol(m).topology().atom(a).iac()+1
	      << "\n\t\tSuggested IACs (atomTypeName; occurence):\n";
	  utils::sort(v, true);
	  
	  for(unsigned int i=0; i<v.size(); ++i){
	    wrn << "\t\t\t" 
		<< setw(4) << v[i].type+1 << " (" 
		<< setw(5) << gffc.atomTypeName(v[i].type) << "; " 
		<< setw(4) << v[i].occurence << ")\n";
	  }
	 
	  inconsistency=true;
	}
	
	// check the mass with the IAC
	ffexp.iac2mass(sys.mol(m).topology().atom(a).iac(), v);
	
	if(!v.size()){
	  if(!inconsistency)
	    wrn << "             Inconsistency with building block file:\n";
	  wrn << "\t\tNo atoms found with IAC "
	      << sys.mol(m).topology().atom(a).iac()+1 << "\n"
	      << "\t\tMaximum IAC in parameter file: " 
	      << gffc.numAtomTypeNames() << "\n";
	  inconsistency=true;
	}
	else{
	  
	  found=false;
	  for(unsigned int i=0; i< v.size(); ++i){
	    if(gffc.findMass(v[i].type) == sys.mol(m).topology().atom(a).mass()) found=true;
	  }
	  if(!found){
	    if(!inconsistency)
	      wrn << "             Inconsistency with building block file:\n";
	    wrn << "\t\tNo atoms found with IAC " 
		<< sys.mol(m).topology().atom(a).iac()+1 << " and Mass "
		<< sys.mol(m).topology().atom(a).mass()
		<< "\n\t\tSuggested Masstype (mass; occurence):\n";
	    utils::sort(v, true);
	    
	    for(unsigned int i=0; i<v.size(); ++i){
	      wrn << "\t\t\t"
		  << setw(4) << v[i].type << " (" 
		  << setw(9) << gffc.findMass(v[i].type) << "; " << setw(5) 
		  << v[i].occurence << ")";
	    }
	    inconsistency=true;
	  }
	  
	  // check the charge with the IAC
	  ffexp.iac2charge(sys.mol(m).topology().atom(a).iac(), v);
	  
	  found=false;
	  for(unsigned int i=0; i< v.size(); ++i){
	    if(ffexp.charge(v[i].type) == 
	       sys.mol(m).topology().atom(a).charge()) found=true;
	  }
	  if(!found){
	    if(!inconsistency)
	      wrn << "             Inconsistency with building block file:\n";
	    wrn << "\t\tNo atoms found with IAC " 
		<< sys.mol(m).topology().atom(a).iac()+1 << " and charge "
		<< sys.mol(m).topology().atom(a).charge()
		<< "\n\t\tSuggested Charge (occurence):\n";

	    utils::sort(v, false);	    
	    for(unsigned int i=0; i<v.size(); ++i){
	      wrn << "\t\t\t"
		  << setw(9) << ffexp.charge(v[i].type) << " (" 
		  << setw(5) << v[i].occurence << ")\n";
	    }
	    inconsistency=true;
	  }
	}
	if(!inconsistency) wrn << "             OK";
	wrn << endl;
        cout << wrn.str();
	
      }
     
      // now do the bonds
      cout << endl << numbonds[m] << " BONDS :" << endl << endl;
      cout << setw(4) << "mol"
	   << setw(10) << "atom-"
	   << setw(12) << "atom-" 
	   << setw(12) << "atom-"
	   << setw(8) << "bond-" 
	   << setw(6) << "found" 
	   << setw(8) << "most"
	   << "  alternatives"<< endl;
      
      cout << setw(4) << "# "
	   << setw(10) << "numbers"
	   << setw(12) << "names"
	   << setw(12) << "IAC"
	   << setw(8) << "type"
	   << setw(14) << "common"
	   << endl;   
      BondIterator bi(sys.mol(m).topology());
      
      for(;bi;++bi){
	int type=bi().type();
	cout << setw(4) << m+1;
	cout << setw(5) << bi()[0]+1 << "-" << setw(4) << bi()[1]+1
	     << setw(7) << sys.mol(m).topology().atom(bi()[0]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(bi()[1]).name()
	     << setw(7) << sys.mol(m).topology().atom(bi()[0]).iac()+1 << "-"
	     << setw(4) << sys.mol(m).topology().atom(bi()[1]).iac()+1
	     << setw(8) << type+1;
	


	// create bond of the types
	Bond b(sys.mol(m).topology().atom(bi()[0]).iac(),
	       sys.mol(m).topology().atom(bi()[1]).iac());
	ffexp.iac2bond(b,v);

	bool found=false;
	bool first=false;
	utils::sort(v, false);
	
	for(unsigned int i=0; i< v.size(); ++i){
	  if(v[i].type == type) {
	    found=true;
	    if(i==0) first=true;
	  }
	}

	int nalt=v.size();
	if(found) {
	  cout << setw(6) << "X";
	  --nalt;
	}
	else cout << setw(6) << " ";
	if(first) cout << setw(6) << "X";
	else cout << setw(6) << " ";
	cout << setw(6) << nalt << endl;
	if(nalt){
	  cout << "\t\t" << setw(4) << "type "
	       << setw(16) << "force constant"
	       << setw(14) << "bond length"
	       << setw(16) << "(occurrence)\n";
	
	  for(unsigned int i=0; i<v.size(); ++i){
	    cout << "\t\t" 
		 << setw(4) << v[i].type+1 << ": ";
	    cout.precision(3);
	    cout.setf(ios::scientific, ios::floatfield);
	    
	    cout << setw(15) << gffc.bondType(v[i].type).fc();
	    cout.setf(ios::fixed, ios::floatfield);
	    cout << setw(14) << gffc.bondType(v[i].type).b0();
	    cout << "   (" 
		 << setw(5) << v[i].occurence << ")\n";
	  }
	}
	
	
	
      }      
      // now do the angles
      cout << endl << numangles[m] << " ANGLES :" << endl << endl;
      cout << setw(4) << "mol"
	   << setw(15) << "atom-"
	   << setw(17) << "atom-" 
	   << setw(17) << "atom-"
	   << setw(8) << "angle-" 
	   << setw(6) << "found" 
	   << setw(8) << "most"
	   << "  alternatives"<< endl;
      cout << setw(4) << "# "
	   << setw(15) << "numbers"
	   << setw(17) << "names"
	   << setw(17) << "IAC"
	   << setw(8) << "type"
	   << setw(14) << "common"
	   << endl;   
      AngleIterator ai(sys.mol(m).topology());
      
      for(;ai;++ai){
	int type=ai().type();
	cout << setw(4) << m+1;
	cout << setw(5) << ai()[0]+1 << "-" << setw(4) << ai()[1]+1
	     << "-" << setw(4) << ai()[2]+1
	     << setw(7) << sys.mol(m).topology().atom(ai()[0]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(ai()[1]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(ai()[2]).name()
	     << setw(7) << sys.mol(m).topology().atom(ai()[0]).iac()+1 << "-"
	     << setw(4) << sys.mol(m).topology().atom(ai()[1]).iac()+1 << "-"
	     << setw(4) << sys.mol(m).topology().atom(ai()[2]).iac()+1
	     << setw(8) << type+1;

	// create angle of the types
	Angle a(sys.mol(m).topology().atom(ai()[0]).iac(),
		sys.mol(m).topology().atom(ai()[1]).iac(),
		sys.mol(m).topology().atom(ai()[2]).iac());
	ffexp.iac2angle(a,v);

	bool found=false;
        bool first=false;
	utils::sort(v, false);
	
	for(unsigned int i=0; i< v.size(); ++i){
	  if(v[i].type == type) {
	    found=true;
	    if(i==0) first=true;
	  }
	}

	int nalt=v.size();
	if(found) {
	  cout << setw(6) << "X";
	  --nalt;
	}
	else cout << setw(6) << " ";
	if(first) cout << setw(6) << "X";
	else cout << setw(6) << " ";
	cout << setw(6) << nalt << endl;
	if(nalt){
	  cout << "\t\t" << setw(4) << "type "
	       << setw(16) << "force constant"
	       << setw(14) << "angle"
	       << setw(16) << "(occurrence)\n";
	  
	  for(unsigned int i=0; i<v.size(); ++i){
	    cout << "\t\t" 
		 << setw(4) << v[i].type+1 << ": ";
	    cout.precision(3);
	    cout.setf(ios::scientific, ios::floatfield);
	    
	    cout << setw(15) << gffc.angleType(v[i].type).fc();
	    cout.setf(ios::fixed, ios::floatfield);
	    cout << setw(14) << gffc.angleType(v[i].type).t0();
	    cout << "   (" 
		 << setw(5) << v[i].occurence << ")\n";
	  }
	}
      }
      
      // now do the impropers
      cout << endl << numimp[m] << " IMPROPER DIHEDRALS :" << endl << endl;
      cout << setw(4) << "mol"
	   << setw(20) << "atom-"
	   << setw(22) << "atom-" 
	   << setw(22) << "atom-"
	   << setw(12) << "improper-" 
	   << setw(6) << "found" 
	   << setw(8) << "most"
	   << "  alternatives"<< endl;
      cout << setw(4) << "# "
	   << setw(20) << "numbers"
	   << setw(22) << "names"
	   << setw(22) << "IAC"
	   << setw(12) << "type"
	   << setw(14) << "common"
	   << endl;   
      ImproperIterator ii(sys.mol(m).topology());
      
      for(;ii;++ii){
	int type=ii().type();
	cout << setw(4) << m+1;
	cout << setw(5) << ii()[0]+1 << "-" << setw(4) << ii()[1]+1 << "-"
	     << setw(4) << ii()[2]+1 << "-" << setw(4) << ii()[3]+1
	     << setw(7) << sys.mol(m).topology().atom(ii()[0]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(ii()[1]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(ii()[2]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(ii()[3]).name()
	     << setw(7) << sys.mol(m).topology().atom(ii()[0]).iac()+1 << "-"
	     << setw(4) << sys.mol(m).topology().atom(ii()[1]).iac()+1 << "-"
	     << setw(4) << sys.mol(m).topology().atom(ii()[2]).iac()+1 << "-"
	     << setw(4) << sys.mol(m).topology().atom(ii()[3]).iac()+1
	     << setw(12) << type+1;

	// create angle of the types
	Improper i(sys.mol(m).topology().atom(ii()[0]).iac(),
		   sys.mol(m).topology().atom(ii()[1]).iac(),
		   sys.mol(m).topology().atom(ii()[2]).iac(),
		   sys.mol(m).topology().atom(ii()[3]).iac());
	ffexp.iac2improper(i,v);

	bool found=false;
        bool first=false;
	utils::sort(v, false);
	
	for(unsigned int i=0; i< v.size(); ++i){
	  if(v[i].type == type) {
	    found=true;
	    if(i==0) first=true;
	  }
	}
	int nalt=v.size();
	if(found) {
	  cout << setw(6) << "X";
	  --nalt;
	}
	else cout << setw(6) << " ";
	if(first) cout << setw(6) << "X";
	else cout << setw(6) << " ";
	cout << setw(6) << nalt << endl;
	if(nalt){
	  cout << "\t\t" << setw(4) << "type "
	       << setw(16) << "force constant"
	       << setw(14) << "improper"
	       << setw(16) << "(occurrence)\n";
	  
	  for(unsigned int i=0; i<v.size(); ++i){
	    cout << "\t\t" 
		 << setw(4) << v[i].type+1 << ": ";
	    cout.precision(3);
	    cout.setf(ios::scientific, ios::floatfield);
	    
	    cout << setw(15) << gffc.improperType(v[i].type).fc();
	    cout.setf(ios::fixed, ios::floatfield);
	    cout << setw(14) << gffc.improperType(v[i].type).q0();
	    cout << "   (" 
		 << setw(5) << v[i].occurence << ")\n";
	  }
	}
      }
      // now do the dihedrals
      cout << endl << numdih[m] << " DIHEDRAL ANGLES :" << endl << endl;
      cout << setw(4) << "mol"
	   << setw(20) << "atom-"
	   << setw(22) << "atom-" 
	   << setw(22) << "atom-"
	   << setw(12) << "dihedral-" 
	   << setw(6) << "found" 
	   << setw(8) << "most"
	   << "  alternatives"<< endl;
      cout << setw(4) << "# "
	   << setw(20) << "numbers"
	   << setw(22) << "names"
	   << setw(22) << "IAC"
	   << setw(12) << "type"
	   << setw(14) << "common"
	   << endl;   
      DihedralIterator di(sys.mol(m).topology());
      
      for(;di;++di){
	int type=di().type();
	cout << setw(4) << m+1;
	cout << setw(5) << di()[0]+1 << "-" << setw(4) << di()[1]+1 << "-"
	     << setw(4) << di()[2]+1 << "-" << setw(4) << di()[3]+1
	     << setw(7) << sys.mol(m).topology().atom(di()[0]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(di()[1]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(di()[2]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(di()[3]).name()
	     << setw(7) << sys.mol(m).topology().atom(di()[0]).iac()+1 << "-"
	     << setw(4) << sys.mol(m).topology().atom(di()[1]).iac()+1 << "-"
	     << setw(4) << sys.mol(m).topology().atom(di()[2]).iac()+1 << "-"
	     << setw(4) << sys.mol(m).topology().atom(di()[3]).iac()+1
	     << setw(12) << type+1;

	// create angle of the types
	Dihedral i(sys.mol(m).topology().atom(di()[0]).iac(),
		   sys.mol(m).topology().atom(di()[1]).iac(),
		   sys.mol(m).topology().atom(di()[2]).iac(),
		   sys.mol(m).topology().atom(di()[3]).iac());
	ffexp.iac2dihedral(i,v);

	bool found=false;
        bool first=false;
	utils::sort(v, false);
	
	for(unsigned int i=0; i< v.size(); ++i){
	  if(v[i].type == type) {
	    found=true;
	    if(i==0) first=true;
	  }
	}
	int nalt=v.size();
	if(found) {
	  cout << setw(6) << "X";
	  --nalt;
	}
	else cout << setw(6) << " ";
	if(first) cout << setw(6) << "X";
	else cout << setw(6) << " ";
	cout << setw(6) << nalt << endl;
	if(nalt){
	  cout << "\t\t" << setw(4) << "type "
	       << setw(16) << "force constant"
	       << setw(14) << "phase shift"
	       << setw(14) << "multiplicity"
	       << setw(16) << "(occurrence)\n";
	    
	  for(unsigned int i=0; i<v.size(); ++i){
	    cout << "\t\t" 
		 << setw(4) << v[i].type+1 << ": ";
	    cout.precision(3);
	    cout.setf(ios::scientific, ios::floatfield);
	    
	    cout << setw(15) << gffc.dihedralType(v[i].type).fc();
	    cout.setf(ios::fixed, ios::floatfield);
	    cout.precision(1);
	    cout << setw(14) << gffc.dihedralType(v[i].type).pd();
	    cout << setw(14) << gffc.dihedralType(v[i].type).np();
	    
	    cout << "   (" 
		 << setw(5) << v[i].occurence << ")\n";
	  }
	}
      }
    }
    
  }
  
  if(args.count("coord")>0){
    cout << endl << "Read in coordinates and calculated covalent energies:"
	 << endl;
    
    int count=0;
    int type=0;
    double totbonds[nummol], totangles[nummol], totimp[nummol], totdih[nummol];
    
    // loop over the properties once again to print
    // bonds
    cout << endl << tnumbonds << " BONDS :" << endl << endl;
    cout << setw(4) << "mol"
	 << setw(10) << "atom-"
	 << setw(12) << "atom-" 
	 << setw(13) << "force-"
	 << setw(10) << "b0"
	 << setw(16) << "b in x"
	 << setw(16) << "energy" << endl;
    cout << setw(4) << "# "
	 << setw(10) << "numbers"
	 << setw(12) << "names"
	 << setw(13) << "constant"
	 << endl;  
    for(int m=0; m<sys.numMolecules(); m++){
      BondIterator bi(sys.mol(m).topology());
      totbonds[m]=0;
      
      for(;bi;++bi){
	type=bi().type();
	cout << setw(4) << m+1;
	cout << setw(5) << bi()[0]+1 << "-" << setw(4) << bi()[1]+1
	     << setw(7) << sys.mol(m).topology().atom(bi()[0]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(bi()[1]).name();
	cout.precision(3);
	cout.setf(ios::scientific, ios::floatfield);
	
	cout << setw(13) << gff.bondType(type).fc();
	cout.setf(ios::fixed, ios::floatfield);
	cout << setw(10) << gff.bondType(type).b0();
	cout.precision(5);
	cout << setw(16) << props[count]->getValue();
	cout.setf(ios::scientific, ios::floatfield);
	cout << setw(16) << en.cov(count) << endl;
	totbonds[m]+=en.cov(count);
	
	count++;
	
      }
    }
    // angles
    cout << endl << tnumangles << " ANGLES :" << endl << endl;
    cout << setw(4) << "mol"
	 << setw(15) << "atom-"
	 << setw(17) << "atom-" 
	 << setw(13) << "force-"
	 << setw(10) << "b0"
	 << setw(16) << "b in x"
	 << setw(16) << "energy" << endl;
    cout << setw(4) << "# "
	 << setw(15) << "numbers"
	 << setw(17) << "names"
	 << setw(13) << "constant"
	 << endl;  
    for(int m=0; m<sys.numMolecules(); m++){
      totangles[m]=0;
      AngleIterator bi(sys.mol(m).topology());
      
      for(;bi;++bi){
	type=bi().type();
	cout << setw(4) << m+1;
	cout << setw(5) << bi()[0]+1 << "-" << setw(4) << bi()[1]+1
	     << "-" << setw(4) << bi()[2]+1
	     << setw(7) << sys.mol(m).topology().atom(bi()[0]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(bi()[1]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(bi()[2]).name();
	cout.precision(3);
	cout.setf(ios::scientific, ios::floatfield);
	
	cout << setw(13) << gff.angleType(type).fc();
	cout.setf(ios::fixed, ios::floatfield);
	cout << setw(10) << gff.angleType(type).t0();
	cout.precision(5);
	cout << setw(16) << props[count]->getValue();
	cout.setf(ios::scientific, ios::floatfield);
	cout << setw(16) << en.cov(count) << endl;
	totangles[m]+=en.cov(count);
	
	count++;
	
      }
    }
    // impropers
    cout << endl << tnumimp << " IMPROPER DIHEDRALS :" << endl << endl;
    cout << setw(4) << "mol"
	 << setw(20) << "atom-"
	 << setw(22) << "atom-" 
	 << setw(13) << "force-"
	 << setw(10) << "b0"
	 << setw(16) << "b in x"
	 << setw(16) << "energy" << endl;
    cout << setw(4) << "# "
	 << setw(20) << "numbers"
	 << setw(22) << "names"
	 << setw(13) << "constant"
	 << endl;  
    for(int m=0; m<sys.numMolecules(); m++){
      totimp[m]=0;
      ImproperIterator bi(sys.mol(m).topology());
      
      for(;bi;++bi){
	type=bi().type();
	cout << setw(4) << m+1;
	cout << setw(5) << bi()[0]+1 << "-" << setw(4) << bi()[1]+1 << "-"
	     << setw(4) << bi()[2]+1 << "-" << setw(4) << bi()[3]+1
	     << setw(7) << sys.mol(m).topology().atom(bi()[0]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(bi()[1]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(bi()[2]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(bi()[3]).name();
	cout.precision(3);
	cout.setf(ios::scientific, ios::floatfield);
	
	cout << setw(13) << gff.improperType(type).fc();
	cout.setf(ios::fixed, ios::floatfield);
	cout << setw(10) << gff.improperType(type).q0();
	cout.precision(5);
	cout << setw(16) << props[count]->getValue();
	cout.setf(ios::scientific, ios::floatfield);
	cout << setw(16) << en.cov(count) << endl;
	totimp[m]+=en.cov(count);
	
	count++;
	
      }
    }
    // dihedrals
    cout << endl << tnumdih << " DIHEDRAL ANGLES :" << endl << endl;
    cout << setw(4) << "mol"
	 << setw(20) << "atom-"
	 << setw(22) << "atom-" 
	 << setw(13) << "force-"
	 << setw(6)  << "pd"
	 << setw(4)  << "np"
	 << setw(16) << "b in x"
	 << setw(16) << "energy" << endl;
    cout << setw(4) << "# "
	 << setw(20) << "numbers"
	 << setw(22) << "names"
	 << setw(13) << "constant"
	 << endl;  
    for(int m=0; m<sys.numMolecules(); m++){
      totdih[m]=0;
      DihedralIterator bi(sys.mol(m).topology());
      
      for(;bi;++bi){
	type=bi().type();
	cout << setw(4) << m+1;
	cout << setw(5) << bi()[0]+1 << "-" << setw(4) << bi()[1]+1 << "-"
	     << setw(4) << bi()[2]+1 << "-" << setw(4) << bi()[3]+1
	     << setw(7) << sys.mol(m).topology().atom(bi()[0]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(bi()[1]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(bi()[2]).name() << "-"
	     << setw(4) << sys.mol(m).topology().atom(bi()[3]).name();
	cout.precision(3);
	cout.setf(ios::scientific, ios::floatfield);
	
	cout << setw(13) << gff.dihedralType(type).fc();
	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(1);
	cout << setw(6) << gff.dihedralType(type).pd();
	cout << setw(4) << gff.dihedralType(type).np();
	
	cout.precision(5);
	cout << setw(16) << props[count]->getValue();
	cout.setf(ios::scientific, ios::floatfield);
	cout << setw(16) << en.cov(count) << endl;
	totdih[m]+=en.cov(count);
	
	count++;
	
      }
    }
    // now summarize some energies
    cout.setf(ios::scientific, ios::floatfield);
    double ttbonds=0, ttangles=0, ttimp=0, ttdih=0;
    
    cout << endl << "SUMMARY :" << endl << endl;
    cout << "Total energies" << endl;
    cout << setw(10) << "molecule"
	 << setw(15) << "bonds"
	 << setw(15) << "angles"
	 << setw(15) << "impropers"
	 << setw(15) << "dihedrals" 
	 << setw(15) << "total" << endl;
    for(int m=0; m<nummol; m++){
      cout << setw(10) << m+1
	   << setw(15) << totbonds[m]
	   << setw(15) << totangles[m]
	   << setw(15) << totimp[m]
	   << setw(15) << totdih[m] 
	   << setw(15) << totbonds[m]+totangles[m]+totimp[m]+totdih[m]<< endl;
      ttbonds+=totbonds[m];
      ttangles+=totangles[m];
      ttimp+=totimp[m];
      ttdih+=totdih[m];
    }
    cout << endl;
    
    if(nummol>1)
      cout << setw(10) << "total"
	   << setw(15) << ttbonds
	   << setw(15) << ttangles
	   << setw(15) << ttimp
	   << setw(15) << ttdih 
	   << setw(15) << ttbonds+ttangles+ttimp+ttdih << endl;
    cout << setw(10) << "average"
	 << setw(15) << ttbonds/tnumbonds
	 << setw(15) << ttangles/tnumangles
	 << setw(15) << ttimp/tnumimp
	 << setw(15) << ttdih/tnumdih << endl;
    cout << endl << "Total covalent energy: " << en.cov() << endl;

    if(error)
      cout << endl << "There were " << error << " warnings about the charge"
	   << " groups, exclusions, bonds, angles or improper dihedrals" 
	   << endl;
    
  }
}
 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}







