/* pdb2g96.cc  This program reads in a topology and a pdb file.
 *             it will then try to generate a gromos-coordinate file
 *             with the atoms in the correct order
 *
 * 
 * Only parses the ATOM and HETATM records from the coordinate section.
 * See 
 * http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
 * for a specification of the PDB file format.
 * For ATOMs and HETATMs we find:
*/

#define ATOM substr(0,4)
#define HETATM substr(0,6)
// the following two are not strictly standard
#define ATOMNAME substr(12, 5)
#define RESNAME substr(17, 4)
//
#define RESNUM substr(22, 4)
#define COORDX substr(30, 8)
#define COORDY substr(38, 8)
#define COORDZ substr(46, 8)

#include "../src/args/Arguments.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"

#include <string>
#include <list>
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace gcore;
using namespace gio;
using namespace gmath;

using namespace args;
using namespace std;

/* ugly but functional c++ hack to strip whitespace
 * from a string */
string stripWhite(string s){

  istringstream bla(s.c_str());
  string fasel; 
  bla >> fasel;
  return fasel;
}

/*
 * Reads the ATOM and HETATM lines from a pdb file into
 * a list<vector<string>>, one vector<string> per residue.
*/
list< vector<string> > readPdbAtoms(ifstream &pdbFile){
    
    int resNum = -1;
    string inPdbLine;
    vector<string> pdbResidue;
    list< vector<string> > pdbResidues;

    while(!pdbFile.eof()){
      getline(pdbFile, inPdbLine);
      if(inPdbLine.ATOM == "ATOM" || 
        inPdbLine.HETATM == "HETATM"){

        // check if we're in a new residue
        if(atoi(inPdbLine.RESNUM.c_str()) != resNum){

          resNum = atoi(inPdbLine.RESNUM.c_str());

          // if we're not in the first residue
          if(pdbResidue.size())
            pdbResidues.push_back(pdbResidue);

          pdbResidue.clear();
        }
        pdbResidue.push_back(inPdbLine);
      }
    }
   
    // push the last residue
    pdbResidues.push_back(pdbResidue);

    return pdbResidues;
}

vector<string> nextPdbResidue(list< vector<string> > &pdbResidues){

  vector<string> pdbResidue;

  if(pdbResidues.begin() != pdbResidues.end()) 
    pdbResidue = *(pdbResidues.begin());

  return pdbResidue;
}

void checkResidueName(vector<string> pdbResidue, string resName){

  if(!pdbResidue.size()){
    ostringstream os;
    os << "Error: Emtpy Residue.\n"
       << "No coordinates in pdb file.";
    
    throw gromos::Exception("pdb2g96", os.str());
  }

  if(stripWhite(pdbResidue[0].RESNAME) != resName){
    ostringstream os;
    os << "Error: Residue names do not match.\n"
       << "\tIn topology: " << resName
       << ", in pdb file: " << pdbResidue[0].RESNAME;
    
    throw gromos::Exception("pdb2g96", os.str());
  }
}

void warnNotFoundAtom(int atomNum, string atomName, 
  int resNum, string resName){

  cerr << "Warning: Could not find atom "
       << atomNum + 1
       << " ("
       << atomName
       << ")"
       << ","

       << " in residue "
       << resNum + 1
       << " ("
       << resName
       << ")"
       << ".\n"

       << "\tSet coordinates to (0.0 0.0 0.0)\n";
}

void warnIgnoredAtoms(vector<string> pdbAtoms){

  for(unsigned int lineNum = 0;
    lineNum < pdbAtoms.size();
    lineNum++)

    cerr << "Warning: Ignored atom "
         << stripWhite(pdbAtoms[lineNum].ATOMNAME)

         << " in residue "
         << stripWhite(pdbAtoms[lineNum].RESNUM)
         << " ("
         << stripWhite(pdbAtoms[lineNum].RESNAME)
         << ").\n";
}


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pdb", "out"};
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pdb <pdb coordinates>\n";
  usage += "\t@out <resulting g96-file> (optional, defaults to stdout)\n";
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // open and read pdb file
    ifstream pdbFile(args["pdb"].c_str());
    list< vector<string> > pdbResidues = readPdbAtoms(pdbFile);
    // loop over all molecules
    for(int molNum = 0; 
      molNum < sys.numMolecules(); 
      molNum++){

      // loop over all residues
      int firstAtomNum = 0, lastAtomNum = 0;
      int resNum = 0;
      for(int thisResNum = 0; 
        thisResNum < sys.mol(molNum).topology().numRes();
        thisResNum++ , resNum++){

        vector<string> pdbResidue = nextPdbResidue(pdbResidues);
        // if the residues in the pdb and the topology are
        // not identical, skip this loop.
        try{
          checkResidueName(pdbResidue, 
            sys.mol(molNum).topology().resName(thisResNum));
          pdbResidues.pop_front();
        }
        catch(gromos::Exception &e){
          cerr << e.what() << endl;
          cerr << " Could not read residue number " << resNum + 1;
          cerr << " from pdb file." << endl;
          cerr << "Skipped" << endl;
          continue;
        }
       
        /* 
         * determine the first and the last atom number of 
         * this residue in the topology 
        */
        for(firstAtomNum = lastAtomNum; 
          sys.mol(molNum).topology().resNum(firstAtomNum) != thisResNum;
          firstAtomNum++);

        for(lastAtomNum = firstAtomNum;
          lastAtomNum < sys.mol(molNum).topology().numAtoms() &&
          sys.mol(molNum).topology().resNum(lastAtomNum) == thisResNum; 
          lastAtomNum++);

        /*
         * for every atom in the topology residue, 
         * look for an atom in the pdb residue with the same name, 
         * and import its coordinates. If we can't find one, 
         * set the coords to 0,0,0 and issue a warning.
        */
        for(int atomNum = firstAtomNum; 
          atomNum < lastAtomNum; 
          atomNum++){

          string inPdbLine = "";
          bool foundAtom = false;

          for(unsigned int pdbAtomNum = 0; 
            pdbAtomNum < pdbResidue.size(); 
            pdbAtomNum++){

            inPdbLine = pdbResidue[pdbAtomNum];
  
            if(sys.mol(molNum).topology().atom(atomNum).name() == 
              stripWhite(inPdbLine.ATOMNAME) &&
              !foundAtom){

              foundAtom = true;  
              sys.mol(molNum).pos(atomNum) = Vec(
                0.1 * atof(inPdbLine.COORDX.c_str()),
                0.1 * atof(inPdbLine.COORDY.c_str()),
                0.1 * atof(inPdbLine.COORDZ.c_str())
              );
              pdbResidue.erase(pdbResidue.begin() + pdbAtomNum);
            }
          }
          if(!foundAtom){
	    sys.mol(molNum).pos(atomNum) = Vec(0.0, 0.0, 0.0);
            warnNotFoundAtom(atomNum,
              sys.mol(molNum).topology().atom(atomNum).name(),
              resNum,
              sys.mol(molNum).topology().resName(thisResNum));
	  }
	}
	// print a warning for the pdb atoms that were ignored 
        warnIgnoredAtoms(pdbResidue);
      }
    }
    // This should be it

    // now define an output stream and write the coordinates
    OutG96S oc;
    try{
      args.check("out", 1);
      ofstream fout(args["out"].c_str());
      oc.open(fout);
      oc.writeTitle("pdb2g96: Reordered atoms from " + args["pdb"]);
      oc << sys;
    }
    catch(const gromos::Exception &e){
      oc.open(cout);
      oc.writeTitle("pdb2g96: Reordered atoms from " + args["pdb"]);
      oc << sys;
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
