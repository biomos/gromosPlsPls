//checktopo reads in a topology and a coordinate file and writes out the
//          the energies for all bonded interactions.

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
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Energy.h"
#include "../src/utils/Neighbours.h"
#include "../src/gmath/Vec.h"
#include <vector>
#include <iomanip>
#include <iostream>
#include <strstream>
#include <fstream>

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

  
int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "coord"};
  int nknowns =3;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type> <gather method>\n";
  usage += "\t@coord <coordinate file>\n";
 
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
  double totcharge[nummol];

  // loop over all bonds
  for(int m=0; m<nummol; m++){
    numbonds[m]=0;
    BondIterator bi(sys.mol(m).topology());
    for(;bi;++bi){
      ostrstream os;
      os << "d%" << m+1 << ":" << bi()[0]+1 << "," << bi()[1]+1 << ends;
      props.addSpecifier(os.str());
      numbonds[m]++;
    }
  }
  // loop over all angles
  for(int m=0; m<nummol; m++){
    numangles[m]=0;
    AngleIterator ai(sys.mol(m).topology());
    for(;ai;++ai){
      ostrstream os;
      os << "a%" << m+1 << ":" << ai()[0]+1 << "," << ai()[1]+1 << "," << ai()[2]+1 << ends;
      props.addSpecifier(os.str());
      numangles[m]++;
    }
  }
  // loop over all impropers
  for(int m=0; m<nummol; m++){
    numimp[m]=0;
    ImproperIterator ii(sys.mol(m).topology());
    for(;ii;++ii){
      ostrstream os;
      os << "t%" << m+1 << ":" << ii()[0]+1 << "," << ii()[1]+1 << "," << ii()[2]+1 << "," << ii()[3]+1 << ends;
      props.addSpecifier(os.str());
      numimp[m]++;
    }
  }
  // loop over all dihedrals
  for(int m=0; m<nummol; m++){
    numdih[m]=0;
    DihedralIterator di(sys.mol(m).topology());
    for(;di;++di){
      ostrstream os;
      os << "t%" << m+1 << ":" << di()[0]+1 << "," << di()[1]+1 << "," << di()[2]+1 << "," << di()[3]+1 << ends;
      props.addSpecifier(os.str());
      numdih[m]++;
    }
  }

  // calculate the total charge
  for(int m=0; m<nummol; m++){
    totcharge[m]=0.0;
    for(int a=0; a<sys.mol(m).numAtoms(); a++)
      totcharge[m]+=sys.mol(m).topology().atom(a).charge();
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
       << "Performing some basic checks on the charge groups, bonds,"
       << " angles and improper dihedrals..." << endl;
  int error=0;
  int chargegroup=0;
  double chrg_precision=0.0001;
  double charge=0.0;
  int chargerest;
  
  // loop over the molecules
  for(int m=0; m< sys.numMolecules(); m++){
    for(int a=0; a < sys.mol(m).topology().numAtoms(); a++){
      //calculate charge group charge
      charge+=sys.mol(m).topology().atom(a).charge();
      if(sys.mol(m).topology().atom(a).chargeGroup()){
        chargegroup++;
	chargerest=int(charge/chrg_precision)%int(1.0/chrg_precision);
	if(chargerest){
          error++;
          cout << setw(5) << error << ". Non-integer valued charge in charge "
	       << "group " << chargegroup << ".\n"
	       << "       Ends at atom " << m +1 << ":" << a +1 << " : " 
	       << chrg_precision*chargerest
	       << endl;
	}
	charge=0.0;
      }
      
      Neighbours neigh(sys, m, a);
      set<int> ex13;
      set<int> ex14;
      
      //check for multiple bonds make sets for the exclusions
      for(unsigned int i=0; i< neigh.size(); i++){
	//check for multiple bonds
	if(neigh[i]>a){
	  for(unsigned int j=i+1; j< neigh.size(); j++){
	    if(neigh[i]==neigh[j]){
	      error++;
	      cout << setw(5) << error << ". More than one bond connecting "
		   << "atoms " << a+1 << " and " << neigh[i]+1 
		   << " in molecule " << m+1 << endl;
	    }
	  }
	}
	// get neighbours of neighbours
        ex13.insert(neigh[i]);
	Neighbours n13(sys, m, neigh[i]);
	for(unsigned int j=0; j<n13.size(); j++){
	  ex13.insert(n13[j]);
	  //and their neighbours
	  Neighbours n14(sys, m, n13[j]);
	  for(unsigned int k=0; k< n14.size(); k++){
	    ex14.insert(n14[k]);
	  }
	}
      }
      // check for the standard exclusion rules
      set<int> ignore14;
      for(int i=0; i< sys.mol(m).topology().atom(a).exclusion().size(); i++){
	int j=sys.mol(m).topology().atom(a).exclusion().atom(i);
	if(j>a && ex13.count(j)==0){
	  error++;
	  cout <<setw(5) << error << ". Atom " << m+1 << ":" << j+1 
	       << " is in the exclusion list of atom " << m+1 
	       << ":" << a+1;
	  if(ex14.count(j)){
	    ignore14.insert(j);
	    cout << " but is a 1,4-neighbour. \n" 
		 << "       This is correct if "
		 <<  m+1 << ":" << a+1 << " and " << m+1 << ":" << j+1
		 << " are part of an aromatic system" << endl;
	  }
	  else
	    cout << " but is not a 1,2 or 1,3 neighbour" << endl;
	}
      }
      for(set<int>::const_iterator iter=ex13.begin(), to=ex13.end(); 
	  iter!=to; ++iter){
	if(*iter>a){
	  int found =0;
	  for(int i=0; i<sys.mol(m).topology().atom(a).exclusion().size();i++)
	    if(*iter==sys.mol(m).topology().atom(a).exclusion().atom(i)) 
	      found=1;
	  if(!found){
	    error++;
	    cout << setw(5) << error << ". Atom " << m+1 << ":" << *iter + 1 
		 << " is a 1,2 or 1,3 neighbour of atom "
		 << m+1 << ":" << a+1 << " but is not in its exclusion list" 
		 << endl;
	  }
	}
      }
        // check for 1,4 exclusions
      for(int i=0; i< sys.mol(m).topology().atom(a).exclusion14().size();i++){
	int j=sys.mol(m).topology().atom(a).exclusion14().atom(i);
	if(j>a && ex14.count(j)==0){
	  error++;
	  cout << setw(5) << error << ". Atom " << m+1 << ":" << j+1 
	       << " is in the 1,4 exclusion list of atom " 
	       << m+1 << ":" << a+1 << " but is not a 1,4 neighbour" << endl;
	}
      }
      for(set<int>::const_iterator iter=ex14.begin(), to=ex14.end(); 
	  iter!=to; ++iter){
	if(*iter>a && ex13.count(*iter)==0){
	  int found =0;
	  for(int i=0; i<sys.mol(m).topology().atom(a).exclusion14().size(); 
	      i++)
	    if(*iter==sys.mol(m).topology().atom(a).exclusion14().atom(i)) 
	      found=1;
	  if(!found && ignore14.count(*iter)==0){
	    error++;
	    cout << setw(5) << error << ". Atom " << m+1 << ":" << *iter + 1 
		 << " is a 1,4 neighbour of atom "
		 << m+1 << ":" << a+1 
		 << " but is not in its 1,4 exclusion list" << endl;
	  }
	}
      }  
    
      //check the angles
      if(neigh.size()>1){
	for(unsigned int i=0; i< neigh.size(); i++){
	  for(unsigned int j=i+1; j< neigh.size(); j++){
	    int b,c;
	    if (neigh[i]<neigh[j]){b=neigh[i]; c=neigh[j];}
	    else                  {b=neigh[j]; c=neigh[i];}
	    AngleIterator ai(sys.mol(m).topology());
	    int count=0;
	    for(;ai;++ai)
	      if(ai()[0]==b&&ai()[1]==a&&ai()[2]==c) count++;
	    if(count==0){
	      error++;
	      cout << setw(5) << error << ". No angle in topology for atoms " 
		   << b+1 << "-" << a+1 << "-" << c+1 << " in molecule " 
		   << m+1 << endl;
	    }
	    if(count>1){
	      error++;
	      cout << setw(5) << error << ". More than one angle in topology "
		   << "for atoms " << b+1 << "-" << a+1 << "-" << c+1 
		   << " in molecule " << m+1 << endl;
	    }
	  }
	}
      }
      //check impropers
      if(neigh.size()==3){
	int at[4]={a,neigh[0], neigh[1], neigh[2]};
	ImproperIterator ii(sys.mol(m).topology());
	int count=0;
	for(;ii;++ii){
	  int found[4]={0,0,0,0};
	  for(int i=0;i<4;i++)
	    for(int j=0; j<4; j++)
	      if(ii()[i]==at[j]) found[j]=1;
	  if(found[0]&&found[1]&&found[2]&&found[3]) count++;
	}
	if(count==0){
	  error++;
	  cout << setw(5) << error << ". No improper dihedral in topology " 
	       << "for atoms " << at[0]+1 << "-" << at[1]+1 << "-" << at[2]+1 
	       << "-" << at[3]+1 << " in molecule " << m+1 << endl;
	}
	if(count>1){
	  error++;
	  cout << setw(5) << error << ". More than one dihedral in topology "
	       << "for atoms " << at[0]+1 << "-" << at[1]+1 << "-" << at[2]+1
	       << "-" << at[3]+1 << " in molecule " << m+1 << endl;
	}
      }
    }
  }
  if(error==0)
    cout << "ok" << endl;
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
      cout << endl << "There were " << error << " warnings about the number"
	   << " of bonds, angles or improper dihedrals" << endl;
    
  }
}
 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}







