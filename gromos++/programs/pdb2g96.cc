// pdb2g96.cc  This program reads in a topology and a pdb file.
//             it will then try to generate a gromos-coordinate file
//             with the atoms in the correct order

//future plans for this program
//- read pdb-lines according to pdb standard. I.e. realize that fields in a 
//  pdb-file can be connected. Use getline. (Or more usefull, write an InPdb-
//  class

#include "../src/args/Arguments.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <strstream>

using namespace gcore;
using namespace gio;
using namespace gmath;

using namespace args;
using namespace std;



int main(int argc, char **argv){

  char *knowns[] = {"topo", "pdb", "out"};
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pdb <pdb coordinates>\n";
  usage += "\t@out <resulting g96-file>\n";
  
  


  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // open pdb-file
    ifstream pdb(args["pdb"].c_str());
    string sdum;
    while(sdum!="ATOM"&&sdum!="HETATM") pdb >> sdum;

    // create vectors to store temporary pdb-coordinates
    vector<string> a_name;
    vector<Vec> a_crd;
    int resno, oldresno;
    string r_name;
    
    //read in the first ATOM line
    pdb >> sdum >> sdum;
    if(sdum.length()>3){
      a_name.push_back(sdum.substr(0,4));
      r_name=sdum.substr(4,sdum.size());
    }
    else{
      a_name.push_back(sdum);
      pdb >> r_name;
    }
    pdb >> resno;
    oldresno = resno;
    Vec tmp;
    pdb >> tmp[0] >> tmp[1] >> tmp[2];
    a_crd.push_back(tmp);
    while(sdum!="ATOM"&&sdum!="HETATM") pdb >> sdum;
    
    
    // loop over all molecules
    for(int m=0;m < sys.numMolecules();m++){

      //loop over all residues
      int nr_atoms=0;
      for(int r=0; r< sys.mol(m).topology().numRes();r++){
        string resname;
	
	//find out how many atoms in this residue
	int nr_res=0;
	while(nr_atoms<sys.mol(m).numAtoms()&&
              sys.mol(m).topology().resNum(nr_atoms)==r) {
          nr_atoms++; 
          nr_res++;
        }
        resname = sys.mol(m).topology().resName(
	 	  sys.mol(m).topology().resNum(nr_atoms-nr_res));
	if(resname!=r_name){
	  ostrstream os;
	  os << "Residue names do not match: In topology "
             << r+1 << " is a " << resname << " but in pdb " << resno
             << " is a " << r_name
             << ends;
	    
          throw gromos::Exception("pbc2g96", os.str());
	}
	
        //read in the complete residue from the pdb
        while(resno==oldresno&&!pdb.eof()){
          pdb >> sdum >> sdum;
          if(sdum.length()>3){
            a_name.push_back(sdum.substr(0,4));
            r_name=sdum.substr(4,sdum.size());
          }
          else{
            a_name.push_back(sdum);
            pdb >> r_name;
          }
          oldresno = resno;
          pdb >> resno;
	  Vec tmp;
          pdb >> tmp[0] >> tmp[1] >> tmp[2];
          a_crd.push_back(tmp);
	  
          while(sdum!="ATOM"&&sdum!="HETATM"&&!pdb.eof()) pdb >> sdum;
	}
	
	// the last entry in the vectors is now the first of the next residue

        
	// we can now follow the topology and sort the pdb-coords
        for(int i=nr_atoms-nr_res;i<nr_atoms;i++){
          string atomname=sys.mol(m).topology().atom(i).name();
	  vector<string>::iterator nm=a_name.begin(), tnm=a_name.end();
	  vector<Vec>::iterator vc=a_crd.begin();
	  if(!pdb.eof()) tnm--;
	  while(nm<tnm&& *nm !=atomname){
	    nm++;
	    vc++;
	  }
          if(*nm==atomname){
	    //we found the atom
            sys.mol(m).pos(i) = 0.1*(*vc);
	    a_name.erase(nm);
	    a_crd.erase(vc);
	  }
	  else{
	    //write a warning that we did not find the atom
            cout << "Could not find atom " << i+1 << " (" << atomname
                 << ") in residue " << r+1 << " (" << resname 
                 << "), molecule " << m+1 << endl;
	    cout << "\tSet coordinates to (0.0 0.0 0.0)\n";
	    sys.mol(m).pos(i) = Vec(0.0,0.0,0.0);
	  }
	}
	int sub=1, sz=a_name.size();
        if(pdb.eof()) sub=0;
	//now, we remove all but the last in our vectors
        for(int q=0;q<sz-sub;q++){
          cout << "Ignored atom " << a_name[0] << " in residue "
               << oldresno << " (" << resname << ") from pdb, molecule "
               << m+1 << endl;
	  a_name.erase(a_name.begin());
	  a_crd.erase(a_crd.begin());
	}
	oldresno = resno;
      }
    }
    // This should be it

    // now define an output stream and write the coordinates
    OutG96S oc;
    ofstream fout(args["out"].c_str());
    
    oc.open(fout);
    oc.writeTitle("pdb2g96: reordered atoms from "+args["pdb"]);
    oc << sys;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

