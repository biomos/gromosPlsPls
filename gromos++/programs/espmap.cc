//espmap calculates the vacuum
//electrostatic potential around a specified
//range of molecules in the topology

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/Box.h"
#include "../src/gio/OutPdb.h"

#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;

int main(int argc, char **argv){
  
  char *knowns[] = {"topo", "pbc", "mol", "grspace", "traj"};
  int nknowns = 5;
  
  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@mol <number of molecules to consider>\n";
  usage += "\t@grspace <grid spacing (default: 0.2 nm)\n";
  usage += "\t@traj <trajectory files>\n";
  
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    
    //  read topology
    InTopology it(args["topo"]);
    //make system
    System sys(it.system());
    
    // which molecules considered?
    vector<int> mols;
    if(args.lower_bound("mol")==args.upper_bound("mol"))
      for(int i=0;i<sys.numMolecules();++i)
	mols.push_back(i);
    else
      for(Arguments::const_iterator it=args.lower_bound("mol");
	  it!=args.upper_bound("mol");++it){
	if(atoi(it->second.c_str())>sys.numMolecules())
	  throw Arguments::Exception(usage);
	mols.push_back(atoi(it->second.c_str())-1);
      }
    
    // get grid spacing
    double space=0.2;
    {
      Arguments::const_iterator iter=args.lower_bound("grspace");
      if(iter!=args.upper_bound("grspace")){
	space=atof(iter->second.c_str());
      }
    }
    
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
    
    // define input coordinate
    InG96 ic;
    Vec center (0.0, 0.0, 0.0), cog (0.0,0.0,0.0), bx(0.0,0.0,0.0),
      grmin(0.0,0.0,0.0), grmax(0.0,0.0,0.0);
    vector<Vec> espgrid;
    vector<double> esp;
    double esptmp=0;
    int atoms=0, nx=0,ny=0,nz=0, numFrames=0;
    double ONE_4PI_EPS0 = 138.9354;
    
    OutCoordinates *oc;
    oc = new OutPdb();
    string ext = ".pdb";
    string extpl = ".pl";
    
    
    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){
      
      // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()){
	numFrames+=1;
	ic >> sys;
	(*pbc.*gathmethod)();
	
	// calculate the center of geometry of the molecules
	
	for(int j=0;j < int (mols.size());++j)
	  for (int i=0;i < sys.mol(mols[j]).numAtoms(); i++) {
	    cog = cog + sys.mol(mols[j]).pos(i);
	    ++atoms;
	  }
	
	cog = (1.0/double(atoms))*cog;
	
	
	//build grid positions
	Box box = sys.box();
	bx[0]= rint(box[0]);bx[1]= rint(box[1]);bx[2]= rint(box[2]);
	
	nx=int(rint(bx[0]/space)); ny=int (rint(bx[1]/space)); nz=int (rint(bx[2]/space));
	//make sure we have equal number of points in x,y,z
	//dont know if this is nescessary for the programs that read in grids
	//i guess it wont hurt either...
	nx=max(nx,ny); nx=max(nx,nz);
	ny=nx; nz=ny;
	
	Vec start = cog-bx/2;
	
	
	for(int i=0;i<nx;i++){
	  for(int j=0;j<ny;j++){
	    for(int k=0;k<nz;k++){ 
	      Vec tmp;      
	      tmp[0]=start[0]+space*i;
	      tmp[1]=start[1]+space*j;
	      tmp[2]=start[2]+space*k;
	      espgrid.push_back(tmp);
	      
	    }
	  }
	}
	grmin=start*10;
	grmax=espgrid[espgrid.size()-1]*10;
	
	
	//calculate esp from gromos charges
	
	for (int i=0;i < int (espgrid.size()); ++i){
	  Vec tmp = espgrid[i];
	  for (int j=0; j<int (mols.size());++j){
	    for (int i=0;i < sys.mol(mols[j]).numAtoms(); i++) { 
	      Vec rvec = (sys.mol(mols[j]).pos(i))-tmp;
	      double r = rvec.abs();
	      esptmp+=ONE_4PI_EPS0*(sys.mol(mols[j]).topology().atom(i).charge())/r;
	    }
	  }
	  esp.push_back(esptmp);
	  esptmp=0;
	}
	
	
	char outFile[]="FRAME";
	ostringstream out;
	ostringstream outpl;
	if (numFrames < 10){
	  out << outFile <<"_"<<"0000"<< numFrames<<ext;
	  outpl << outFile <<"_"<<"0000"<< numFrames<<extpl;
	}
	else if (numFrames < 100){
	  out << outFile<<"_"<<"000"<< numFrames<<ext;
	  outpl << outFile<<"_"<<"000"<< numFrames<<extpl;
	}
	else if (numFrames < 1000){
	  out << outFile <<"_"<<"00"<< numFrames<<ext;
	  outpl << outFile <<"_"<<"00"<< numFrames<<extpl;
	}
	else {
	  out << outFile <<"_"<<"0"<< numFrames<<ext;
	  outpl << outFile <<"_"<<"0"<< numFrames<<extpl;
	}
	//write pdb
	ofstream os(out.str().c_str()); 
	oc->open(os);       
	oc-> select("SOLUTE");
	oc->writeTitle(out.str());
	
	*oc << sys;
	os.close();
	
	//write pl-file
	ofstream opl;
	opl.open(outpl.str().c_str());
	
	opl << "3" << ' ' << "3" << endl;
	opl << nx << ' ' << ny << ' ' << nz << endl;
	//               opl.setf(scientific, floatfield);
	opl << grmin[2] << ' ' << grmax[2] << ' ' 
	    << grmin[1] << ' ' << grmax[1] << ' '
	    << grmin[0] << ' ' << grmax[0] << endl;
	
	for (int i=0,j=0;i < int (esp.size()); ++i, ++j){
	  if (j==2) {opl << endl; j=0;}
	  opl << ' ' << esp[i] << ' ' ; 
	}
	
	
	
	
	opl.close(); 
	
      }
      ic.close();
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

