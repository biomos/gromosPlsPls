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

#include <string>
#include <list>
#include <vector>
#include <map>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gmath/Vec.h"


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
 * checks if two names are the same or should be 
 * considered the same
 */
bool checkName(multimap<string,string> lib, string nameA, string nameB)
{
  nameA=stripWhite(nameA);
  if(nameA==nameB) return true;
  for(multimap<string,string>::const_iterator iter
	  =lib.lower_bound(nameA), to=lib.upper_bound(nameA);
      iter!=to; ++iter)
    if(iter->second == nameB) return true;
  return false;
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

void checkResidueName(vector<string> pdbResidue, string resName, 
		      multimap<string, string> &libRes){

  if(!pdbResidue.size()){
    ostringstream os;
    os << "Error: Emtpy Residue.\n"
       << "No coordinates in pdb file.";
    
    throw gromos::Exception("pdb2g96", os.str());
  }

  if(!checkName(libRes, pdbResidue[0].RESNAME, resName)){
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

void readLibrary(Ginstream &lib, multimap<string, string> &libRes,
		 map<string, multimap<string, string> > &libAtom)
{

  typedef multimap<string,string>::value_type mapType;
  
  std::vector<std::string> buffer;
  std::vector<std::vector<std::string > > content;
  while(!lib.stream().eof()){
    lib.getblock(buffer);
    if(!lib.stream().eof()){
      if(buffer[buffer.size()-1].find("END")!=0)
	throw gromos::Exception("pdb2g96", "Library file " + lib.name() +
				" is corrupted. No END in "+buffer[0]+
				" block. Got\n"
				+ buffer[buffer.size()-1]);

      content.push_back(buffer);
    }    
  }
  // now loop over the content
  std::vector<std::vector<std::string > >::const_iterator 
    iter=content.begin();
  string resNameA, resNameB, atomNameA, atomNameB;
  
  for( ; iter!=content.end(); ++iter){
    if ((*iter)[0]=="RESIDUES"){    
      for(unsigned int i=1; i< (*iter).size()-1; i++){
	std::istringstream linestream((*iter)[i]);
	linestream >> resNameA >> resNameB;
	libRes.insert(mapType(resNameA, resNameB));
      }
      
    }
    else if((*iter)[0]=="ATOMS"){      
      for(unsigned int i=1; i< (*iter).size()-1; i++){
	std::istringstream linestream((*iter)[i]);
	linestream >> resNameA >> atomNameA >> atomNameB;
	libAtom[resNameA].insert(mapType(atomNameA,atomNameB));
      }
    }
    else
      throw gromos::Exception("pdb2g96", 
	     "Don't know how to handle "+(*iter)[0] + "-block in library");
  }
}

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pdb", "out", "lib"};
  int nknowns = 4;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pdb  <pdb coordinates>\n";
  usage += "\t@out  <resulting GROMOS coordinates> (optional, defaults to stdout)\n";
  usage += "\t@lib  <library for atom and residue names>\n";
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // open and read pdb file
    ifstream pdbFile(args["pdb"].c_str());
    list< vector<string> > pdbResidues = readPdbAtoms(pdbFile);

    // read the library file
    std::multimap<std::string, std::string> libRes;
    std::map<std::string, std::multimap<std::string, std::string> > libAtom;
    if(args.count("lib")>0){
      Ginstream lib(args["lib"]);
      readLibrary(lib, libRes, libAtom);
    }
    
    // loop over all molecules
    for(int molNum = 0; 
	molNum < sys.numMolecules(); 
	molNum++){

      // reserve memory for the coordinates
      sys.mol(molNum).initPos();
      
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
            sys.mol(molNum).topology().resName(thisResNum), libRes);
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
	    
            if(checkName(libAtom[sys.mol(molNum).topology().resName(resNum)],
			 inPdbLine.ATOMNAME,
			 sys.mol(molNum).topology().atom(atomNum).name()) 
	       && !foundAtom){

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
    // everything that is left over, should be solvent
    int countSolvent=0;
    
    while (pdbResidues.size()){
      vector<string> pdbResidue = nextPdbResidue(pdbResidues);
      countSolvent++;
      
      try{
	checkResidueName(pdbResidue, "SOLV", libRes);
	pdbResidues.pop_front();
      }
      catch(gromos::Exception &e){
	cerr << e.what() << endl;
	cerr << " Could not read residue number " << countSolvent;
	cerr << " from pdb file." << endl;
	cerr << "Skipped" << endl;
	continue;
      }
      /*
       * for every atom in the topology residue, 
       * look for an atom in the pdb residue with the same name, 
       * and import its coordinates. If we can't find one, 
       * set the coords to 0,0,0 and issue a warning.
       */
      for(int atomNum = 0; 
	  atomNum < sys.sol(0).topology().numAtoms(); 
	  atomNum++){

	string inPdbLine = "";
	bool foundAtom = false;

	for(unsigned int pdbAtomNum = 0; 
	    pdbAtomNum < pdbResidue.size(); 
	    pdbAtomNum++){

	  inPdbLine = pdbResidue[pdbAtomNum];
	    
	  if(checkName(libAtom["SOLV"],
		       inPdbLine.ATOMNAME,
		       sys.sol(0).topology().atom(atomNum).name()) 
	     && !foundAtom){

	    foundAtom = true;  
	    sys.sol(0).addPos(Vec(0.1 * atof(inPdbLine.COORDX.c_str()),
				    0.1 * atof(inPdbLine.COORDY.c_str()),
				    0.1 * atof(inPdbLine.COORDZ.c_str())));
	    pdbResidue.erase(pdbResidue.begin() + pdbAtomNum);
	  }
	}
	if(!foundAtom){
	  sys.sol(0).addPos(Vec(0.0, 0.0, 0.0));
	  warnNotFoundAtom(atomNum,
			   sys.sol(0).topology().atom(atomNum).name(),
			   countSolvent-1, "SOLV");
	}
      }
      // print a warning for the pdb atoms that were ignored 
      warnIgnoredAtoms(pdbResidue);
    }

    //that's really it

    // now define an output stream and write the coordinates
    OutG96S oc;
    try{
      args.check("out", 1);
      ofstream fout(args["out"].c_str());
      oc.open(fout);
      oc.select("ALL");
      oc.writeTitle("pdb2g96: Reordered atoms from " + args["pdb"]);
      oc << sys;
    }
    catch(const gromos::Exception &e){
      oc.open(cout);
      oc.select("ALL");
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
