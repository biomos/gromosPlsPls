//iondens calculates the ion density around
//the first molecule of the trajectory and
//writes out a pdb file with B-factors; sort of derived from espmap...

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/Box.h"
#include "../src/gio/OutPdb.h"

#include <vector>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <functional>


using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace fit;
using namespace std;
using namespace gmath;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "class", "atoms", "mol", "ref", "grspace", "traj", "ions", "thresholds"};
  int nknowns = 10;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@mol <number of molecules to consider>\n";
  usage += "\t@class <classes of atoms to consider>\n";
  usage += "\t@atoms <atoms to consider>\n";
  usage += "\t@ref <reference coordinates>\n";
  usage += "\t@grspace <grid spacing (default: 0.2 nm)\n";
  usage += "\t@ions <ion molecule numbers in topology (default: 1)\n";
  usage += "\t@thresholds <threshold value for occupancy percentages (default: 20 and 5)\n";
  usage += "\t@traj <trajectory files>\n";
  
 
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    //  read topology
    InTopology it(args["topo"]);
    //make system
    System sys(it.system());
    //input coordinate
    InG96 ic;


    // which molecules considered?
    vector<int> mols;
    if(args.lower_bound("mol")==args.upper_bound("mol"))
      for(int i=0;i<sys.numMolecules();++i)
	mols.push_back(i);
    else
      for(Arguments::const_iterator it=args.lower_bound("mol");
	  it!=args.upper_bound("mol");++it){
	if(atoi(it->second.c_str())>sys.numMolecules())
	  throw Arguments::Exception("Molecule number not in topology");
	mols.push_back(atoi(it->second.c_str())-1);
      }

    //which ions to consider
    vector<int> ions;
    if(args.lower_bound("ions")==args.upper_bound("ions"))
      //	for(int i=0;i<sys.numMolecules();++i)
      ions.push_back(1);
    else
      for(Arguments::const_iterator it=args.lower_bound("ions");
	  it!=args.upper_bound("ions");++it){
	if(atoi(it->second.c_str())>sys.numMolecules())
	  throw Arguments::Exception("Ion number not in topology");
	ions.push_back(atoi(it->second.c_str())-1);
      }

    // get grid spacing
    double space=0.2;
    {
      Arguments::const_iterator iter=args.lower_bound("grspace");
      if(iter!=args.upper_bound("grspace")){
	space=atof(iter->second.c_str());
      }
    }

    // get threshold values
    vector<double> thres;
    if(args.lower_bound("thresholds")==args.upper_bound("thresholds")) {
      thres.push_back(20);thres.push_back(2); 
    }
    else
      for(Arguments::const_iterator it=args.lower_bound("thresholds");
	  it!=args.upper_bound("thresholds");++it){
	thres.push_back(atof(it->second.c_str()));
      }

    System refSys(it.system());

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


    
    // parse boundary conditions
    //Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    //Chris: nothing has been read into system, but the pbc are 
    //       based on them; should it here not be refsys that gets gathered?

    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    (*pbc.*gathmethod)();

    delete pbc;

    Reference ref(&refSys);
    ReferenceParser refP(refSys, args, ref);
    refP.add_ref();

    RotationalFit rf(&ref);
    rf.fit(&refSys); 

    OutCoordinates *oc;
    oc = new OutPdb();
    ofstream os("ref.pdb");
    OutPdb opdb(os);
    opdb << refSys;
    os.close();

    //create a system for the average structure, set the coords to 0.0
    System aver(refSys);
    for (int i=0;i < aver.numMolecules(); ++i){
      for (int j=0; j < aver.mol(i).numAtoms(); ++j){
	aver.mol(i).pos(j)=0.0;
      }
    }


    // calculate the center of geometry of the molecules
    Vec center (0.0, 0.0, 0.0), cog (0.0,0.0,0.0), bx(0.0,0.0,0.0),
      grmin(0.0,0.0,0.0), grmax(0.0,0.0,0.0);

    int atoms=0, nx=0,ny=0,nz=0, numFrames=0;

    for(int j=0;j < int (mols.size());++j)
      for (int i=0;i < refSys.mol(mols[j]).numAtoms(); i++) {
	cog = cog + refSys.mol(mols[j]).pos(i);
	++atoms;
      }

    cog = (1.0/double(atoms))*cog;

    //build grid positions
    vector<Vec> densgrid;

    Box box = refSys.box();
    box[0]*=1.2;
    box[1]*=1.2;
    box[2]*=1.2;
    
    nx=int(rint(box[0]/space)); ny=int(rint(box[1]/space)); nz=int(rint(box[2]/space));

    //make sure we have equal number of points in x,y,z
    //dont know if this is nescessary for the programs that read in grids
    //i guess it wont hurt either...
    nx=max(nx,ny); nx=max(nx,nz);
    ny=nx; nz=ny;

    bx[0]=nx*space; bx[1]=ny*space; bx[2]=nz*space;

    Vec start=-bx/2;
 
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
	for(int k=0;k<nz;k++){ 
	  Vec tmp;      
	  tmp[0]=start[0]+space*i;
	  tmp[1]=start[1]+space*j;
	  tmp[2]=start[2]+space*k;
	  densgrid.push_back(tmp);
	}
      }
    }
    grmin=start;
    grmax=densgrid[densgrid.size()-1];

    vector<int> ioncount; for (int i=0; i < int (densgrid.size()); ++i){ioncount.push_back(0);}


    pbc = BoundaryParser::boundary(sys, args);
     
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

	rf.fit(&sys);	   


	//calculate cog again for all selected molecules
	cog[0]=0.0; cog[1]=0.0; cog[2]=0.0;
	atoms=0;
	for(int j=0;j < int (mols.size());++j)
	  for (int i=0;i < sys.mol(mols[j]).numAtoms(); i++) {
	    cog = cog + sys.mol(mols[j]).pos(i);
	    ++atoms;
	  }

	cog = (1.0/double(atoms))*cog;
	
	//calculate nim with respect to cog from above
	for (int i=0; i < int (ions.size()); ++i){
	  sys.mol(ions[i]).pos(0) = pbc->nearestImage(cog,sys.mol(ions[i]).pos(0),sys.box());
	}

	//sum average position
	for (int i=0;i < aver.numMolecules(); ++i){
	  for (int j=0; j < aver.mol(i).numAtoms(); ++j){
	    aver.mol(i).pos(j)+=sys.mol(i).pos(j);
	  }
	}

 
	double rmin = 1000000, r=0; int p=0;
	//determine ion position with respect to grid
	for (int i=0; i< int(ions.size()); ++i){
	  p=0, rmin = 1000000;
	  Vec ion = sys.mol(ions[i]).pos(0);
	  if(fabs(ion[0]-start[0])>bx[0] || 
	     fabs(ion[1]-start[1])>bx[1] ||
	     fabs(ion[2]-start[2])>bx[2])
	    cout << "outside " << ion[0] << " " << ion[1] << " " << ion[2] << endl;
	  
	  for (int j=0; j < int (densgrid.size()); ++j){
	    Vec tmp = densgrid[j];
	    r = (tmp-ion).abs();
	    if (rmin > r) {
	      rmin = r;
	      p = j;}
	  }
	  ioncount[p]+=1;
	}

      }
      ic.close();
    }
    //average and write out average structure
    for (int i=0;i < aver.numMolecules(); ++i){
      for (int j=0; j < aver.mol(i).numAtoms(); ++j){
	aver.mol(i).pos(j)=aver.mol(i).pos(j)/numFrames;
      }
    }

    ofstream oss("aver.pdb");
    OutPdb oopdb(oss);
    oopdb << aver;
    oss.close();


    //get max_element in ioncount
    vector<int>::iterator maxel = std::max_element(ioncount.begin(),ioncount.end());
  
    //normalize * 100
    vector<double> per;
    vector<double> pernf;
    for (int i=0; i < int (ioncount.size()); ++i){ 
      per.push_back(double(double(ioncount[i])/double(*maxel))*100);
      pernf.push_back(double(double(ioncount[i])/double(numFrames))*100);
      //do something stupid...
      if (per[i] == 100) {per[i]=99.99;}
      if (pernf[i] == 100) {pernf[i]=99.99;}
    }

    //write out 'da grid with b-facs
    ofstream perfile; perfile.open("grid.pdb");
    int count=0;
    for (int i=0; i < int (ioncount.size()); ++i){ 
      if (per[i] > thres[0]) {
	++count;
	Vec tmp = densgrid[i];
	perfile.setf(ios::fixed, ios::floatfield);
	perfile.setf(ios::unitbuf);
	perfile.precision(3);

	perfile << "ATOM";
	perfile.setf(ios::right, ios::adjustfield);
	perfile << setw(7) << count;
	perfile.setf(ios::left, ios::adjustfield);
	perfile << "  " <<setw(4) << "AR";
	perfile << setw(4) << "AR";
	perfile.setf(ios::right, ios::adjustfield);
	perfile << setw(5) << "1" << "    "
		<< setw(8) << tmp[0]*10
		<< setw(8) << tmp[1]*10
		<< setw(8) << tmp[2]*10
		<< "  1.00 "; 
	perfile.precision(2);
	perfile << per[i]
		<< endl;
      }
    }
    perfile.close();

    //write out 'da grid with b-facs
    ofstream pernffile; pernffile.open("gridnf.pdb");
    count=0;
    for (int i=0; i < int (ioncount.size()); ++i){ 
      if (pernf[i] > thres[1]) {
	++count;
	Vec tmp = densgrid[i];
	pernffile.setf(ios::fixed, ios::floatfield);
	pernffile.setf(ios::unitbuf);
	pernffile.precision(3);

	pernffile << "ATOM";
	pernffile.setf(ios::right, ios::adjustfield);
	pernffile << setw(7) << count;
	pernffile.setf(ios::left, ios::adjustfield);
	pernffile << "  " <<setw(4) << "AR";
	pernffile << setw(4) << "AR";
	pernffile.setf(ios::right, ios::adjustfield);
	pernffile << setw(5) << "1" << "    "
		  << setw(8) << tmp[0]*10
		  << setw(8) << tmp[1]*10
		  << setw(8) << tmp[2]*10
		  << "  1.00 "; 
	pernffile.precision(2);
	pernffile << pernf[i]
		  << endl;
      }
    }
    pernffile.close();


  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

