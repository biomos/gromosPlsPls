// root mean square fluctuations rmsf.cc

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Vec.h"


#include <string>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <iostream>

using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "class", "pbc", "ref", "mol", "frskip"};
  int nknowns = 7;
  int mol = 0;
  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@mol <molecules to be considered>\n";
  usage += "\t@class <classes of atoms to consider>\n";
  usage += "\t@frskip <total number of frames to be skipped before averaging>\n";
  usage += "\t@ref <reference coordinates>\n";
  usage += "\t@traj <trajectory files>\n";


  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    
    // read reference coordinates...
    InG96 ic;

    try{
      args.check("ref",1);
      ic.open(args["ref"]);
    }
    catch(const Arguments::Exception &){
      args.check("traj",1);
      ic.open(args["traj"]);
    }
    ic >> refSys;
    ic.close();


    // Parse boundary conditions
 Boundary *pbc = BoundaryParser::boundary(refSys, args);
 // parse gather method
  Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // gather reference system
    (*pbc.*gathmethod)();
    delete pbc;
    
    Reference ref(&refSys);

    // Adding references
    int added=0;
    // which molecules considered?
    vector<int> mols;
    if(args.lower_bound("mol")==args.upper_bound("mol"))
      for(int i=0;i<refSys.numMolecules();++i)
	mols.push_back(i);
    else
      for(Arguments::const_iterator it=args.lower_bound("mol");
	  it!=args.upper_bound("mol");++it){
	if(atoi(it->second.c_str())>refSys.numMolecules())
	  throw Arguments::Exception(usage);
	mols.push_back(atoi(it->second.c_str())-1);
      }
  // number of frames to be skipped
   int  frskip=0;
   {
   Arguments::const_iterator iter=args.lower_bound("frskip");
   if(iter!=args.upper_bound("frskip")){
     frskip=atoi(iter->second.c_str());
   }
    }

    // add classes for fit and essential dynamics
    vector<string> name;
    for(Arguments::const_iterator it=args.lower_bound("class");
	it != args.upper_bound("class"); ++it){
    name.push_back(it->second);
      for(vector<int>::const_iterator mol=mols.begin();
	  mol!=mols.end();++mol)
	ref.addClass(*mol,it->second);
      added=1;
    }

    // System for calculation
    System sys(refSys);

    // get atoms from class 
        vector<int> atomlist;
    for (int i=0;i< int (name.size());++i){
      ref.makePosList(sys, mols[0], name[i], atomlist);
    } 
   
     // sort atom numbers
    std::sort(atomlist.begin(), atomlist.end());

    // did we add anything at all?
    if(!added)
      throw Arguments::Exception(usage);

    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);
    RotationalFit rf(&ref);


 int numFrames = 0;
 int size = atomlist.size();
 double avpos[size*3];
 for (int i = 0;i<size*3;++i){avpos[i]=0;}
 double pvec[size*3];
 Matrix cov(size*3,size*3,0.0);

 for(Arguments::const_iterator iter=args.lower_bound("traj");
  iter!=args.upper_bound("traj"); ++iter){
  ic.open(iter->second);
    // loop over all frames
  while(!ic.eof()){
   numFrames++;
   ic >> sys;	
   (*pbc.*gathmethod)();
   rf.fit(&sys);
   if (numFrames >= frskip){
    // calculate average positions, put coords into one array
   for (int i=0, x=0; i<size;++i, x+=2) {
     int j = atomlist[i];
    pvec[i+x]=sys.mol(mol).pos(j)[0];
    pvec[i+1+x]=sys.mol(mol).pos(j)[1];
    pvec[i+2+x]=sys.mol(mol).pos(j)[2];
    avpos[i+x]+=pvec[i+x];
    avpos[i+1+x]+= pvec[i+1+x];
    avpos[i+2+x]+= pvec[i+2+x];
   }  

   //build one triangle of the covariance matrix
   for (int ii=0; ii<size*3;++ii){
    for (int jj=0; jj<=ii;++jj){
     cov(ii,jj)=cov(ii,jj)+(pvec[ii]*pvec[jj]);
    }
   }
   }
  }
  ic.close();
 }
   
  //average positions, write them out
    for (int i=0, x=0, y=0,z=0; i<size*3;++i) {
     y=atomlist[z];
     avpos[i] = avpos[i]/(numFrames-frskip);
     sys.mol(0).pos(y)[x]=avpos[i];
    x+=1; if (x == 3){x=0;z+=1;}
    }


 int NDIM = size*3; 
    for (int ii=0;ii<NDIM;++ii){
     for (int jj=0;jj<=ii;++jj){
      cov(ii,jj)=(cov(ii,jj)/(numFrames-frskip))-(avpos[ii]*avpos[jj]);
     }
    }  
   
  //calculate trace of the symmetric matrix, i.e. the rmsf
    double rmsf=0;    
    for (int z=0, x=0;z<size;++z, x+=2){
      rmsf=((cov(z+x,z+x)+cov(z+x+1,z+x+1)+cov(z+x+2,z+x+2)));
      rmsf=sqrt(rmsf);
      cout << z+1 << ' ' << rmsf << endl;
    }      
       }

  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


