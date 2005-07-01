// rmsdmat.cc

#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/FastRotationalFit.h"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace bound;
using namespace gio;
using namespace utils;
using namespace args;
using namespace fit;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "pbc", "ref", "atomsrmsd", "atomsfit", 
		    "skip", "stride", "human", "precision", "big"}; 

  const int nknowns = 11;

  string usage = argv[0];
  usage += "\n\t@topo     <topology>\n";
  usage += "\t@pbc          <boundary conditions> <gather type>\n";
  usage += "\t@atomsfit     <atomspecifier: atoms to consider for fit\n";
  usage += "\t[@atomsrmsd   <atomspecifier: atoms to consider for rmsd\n";
  usage += "\t                (only specify if different from atomsfit\n";
  usage += "\t[@skip        <skip frames at beginning>]\n";
  usage += "\t[@stride      <use only every step frame>]\n";
  usage += "\t[@human       (write the matrix in human readable form)]\n";
  usage += "\t[@precision   <number of digits in the matrix (default 4)>\n";
  usage += "\t[@big         use if you have more than ~50'000 structures]\n";
  usage += "\t[@ref         <reference coordinates>]\n";
  usage += "\t@traj         <trajectory files>\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    

    // read the fit atoms
    AtomSpecifier fitatoms(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("atomsfit"),
	to=args.upper_bound("atomsfit");
      for( ; iter!=to; ++iter){
	fitatoms.addSpecifier(iter->second);
      }
      if(fitatoms.size()==0)
	throw gromos::Exception("rmsdmat",
				"No fit atoms given\n");
    }

    // read the rmsd atoms
    AtomSpecifier rmsdatoms(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("atomsrmsd"),
	to=args.upper_bound("atomsrmsd");
      for( ;iter!=to; ++iter){
	rmsdatoms.addSpecifier(iter->second);
      }
      if(rmsdatoms.size()==0)
	rmsdatoms = fitatoms;
      
    }
    
    // for which atoms do we want to keep the coordinates
    AtomSpecifier atoms(fitatoms + rmsdatoms);
    
    // if fitatoms != rmsdatoms keep lists of what to do
    vector<bool> fit_spec, rmsd_spec;
    if(atoms.size()!=fitatoms.size() ||atoms.size() != rmsdatoms.size()){
      fit_spec.resize(atoms.size(), true);
      rmsd_spec.resize(atoms.size(), true);
      for(int i=0; i < atoms.size(); ++i){
	if(fitatoms.findAtom(atoms.mol(i), atoms.atom(i))<0)
	  fit_spec[i] = false;
	if(rmsdatoms.findAtom(atoms.mol(i), atoms.atom(i))<0)
	  rmsd_spec[i] = false;
      }
    }

    FastRotationalFit frf(fit_spec, rmsd_spec);
    
    int skip=0;
    if(args.count("skip")>0) skip = atoi(args["skip"].c_str());
    int stride=1;
    if(args.count("stride")>0) stride = atoi(args["stride"].c_str());

    // read the precision
    int precision=10000;
    if(args.count("precision")>0){
      int ii=atoi(args["precision"].c_str());
      precision=1;
      for(int i=0; i<ii; ++i)
	precision*=10;
    }
    
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // create the vector to store the trajectory
    vector< vector < Vec > > traj;
    vector< Vec > frame(atoms.size());
    
    // read reference coordinates...
    InG96 ic;
    if(args.count("ref") > 0){
      ic.open(args["ref"]);
    }
    else{
      ic.open(args.lower_bound("traj")->second);
    }
    ic >> sys;
    ic.close();

    (*pbc.*gathmethod)();

    // calculate the centre of geometry of the relevant atoms
    Vec cog;
    for(int i=0; i < atoms.size(); ++i){
      cog+= *atoms.coord(i);
    }
    cog /= atoms.size();

    // put it in the trajectory
    for(int i=0; i < atoms.size(); ++i){
      frame[i] = *atoms.coord(i) - cog;
    }
    traj.push_back(frame);
    
    int framenum=0;
    
    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
      iter!=args.upper_bound("traj"); ++iter){

      ic.open(iter->second);

      // loop over all frames
      while(!ic.eof()){
	ic >> sys;
	
	if (!((framenum - skip) % stride)){
	  
	  //pbc call
	  (*pbc.*gathmethod)();
	  
	  Vec cog;
	  for(int i=0; i < atoms.size(); ++i){
	    cog+= *atoms.coord(i);
	  }
	  cog /= atoms.size();
	  for(int i=0; i < atoms.size(); ++i){
	    frame[i] = *atoms.coord(i) - cog;
	  }
	  int err = frf.fit(traj[0], frame);
	  if(err) {
	    ostringstream os;
	    os << "Error while fitting to the reference structure\n"
	       << "Error code " << err << " in frame number " << framenum+1;
	    
	    throw gromos::Exception("rmsdmat", os.str());
	  }
	  
          // store coordinates from sys in traj
	  traj.push_back(frame);
        }
	framenum++;
      }
      ic.close();
    }

    // everything is in the thing; create the thingy
    // open a file
    ofstream fout;
    bool human=false;
    if(args.count("human") >= 0){
      fout.open("RMSDMAT.dat");
      human=true;
    }
    else{
      fout.open("RMSDMAT.bin", ios::out | ios::binary);
    }
    
    // make a double loop
    int num = traj.size();
    Matrix rot(3,3,0);
    double rmsd = 0;
    Matrix unit(3,3,0);
    for(size_t i=0; i<3; i++) unit(i,i)=1;

    cout << "Read " << num << " out of " << framenum 
	 << " structures from trajectory" << endl;
    
    if(human){
      fout << "TITLE\n"
	   << "\trmsd-matrix for " << num-1 << " + 1 (ref) = " 
	   << num << " structures\n"
	   << "END\n"
	   << "RMSDMAT\n"
	   << "# number of frames   skip   stride\n"
	   << num << "\t" << skip << "\t" << stride << "\n"
	   << "# precision\n"
	   << precision << "\n";
      

      for(int i=0; i< num; ++i){
	for(int j=i+1; j < num; ++j){
	  if(frf.fit(rot, traj[i], traj[j])){
	    cout << "error rotational fit on frames " << i + 1 << " and " 
		 << j + 1 
		 << "\nfitting to reference structure instead" << endl;
	    rot=unit;
	  }

	  rmsd = frf.rmsd(rot, traj[i], traj[j]);
	  fout << setw(8) << i 
	       << setw(8) << j 
	       << setw(8) << unsigned(rmsd * precision)
	       << endl;
	}
      }
      fout << "END\n";
    }
    else{
      fout.write((char*)&num, sizeof(int));
      fout.write((char*)&skip, sizeof(int));
      fout.write((char*)&stride, sizeof(int));
      fout.write((char*)&precision, sizeof(int));
      
      if(precision < 1e5 && args.count("big") == -1){
	std::cout << "using 'unsigned short' as format" << std::endl;

	typedef unsigned short ushort;
	
	ushort irmsd;
	for(int i=0; i< num; ++i){
	  for(int j=i+1; j < num; ++j){
	    if(frf.fit(rot, traj[i], traj[j])){
	      cout << "error rotational fit on frames " << i + 1 << " and " 
		 << j + 1 
		 << "\nfitting to reference structure instead" << endl;
	      rot=unit;
	    }
	    
	    rmsd = frf.rmsd(rot, traj[i], traj[j]);
	    irmsd=ushort(rmsd*precision);
	    
	    fout.write((char*)&irmsd, sizeof(ushort));
	  }
	}
      }
      else{
	std::cout << "using 'unsigned int' as format" << std::endl;
	unsigned irmsd;
	for(int i=0; i< num; ++i){
	  for(int j=i+1; j < num; ++j){
	    if(frf.fit(rot, traj[i], traj[j])){
	      cout << "error rotational fit on frames " << i + 1 << " and " 
		 << j + 1 
		 << "\nfitting to reference structure instead" << endl;
	      rot=unit;
	    }
	    
	    rmsd = frf.rmsd(rot, traj[i], traj[j]);
	    irmsd=unsigned(rmsd*precision);
	    
	    fout.write((char*)&irmsd, sizeof(unsigned));
	  }
	}
      }	
    }
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
