// cluster.cc

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/OutPdb.h"
#include "../src/gmath/Vec.h"
#include "../src/gio/TrajArray.h"

#include <vector>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

int main(int argc, char **argv){

  const int MAXTRAJSIZ=8002;

  char *knowns[] = {"topo", "traj", "class", "atoms", "pbc", "ref",
                    "mol", "time","mat","step","ctraj","normsd"};
  int nknowns = 12;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@mol <molecules to be considered>\n";
  usage += "\t@class <classes of atoms to consider>\n";
  usage += "\t@atoms <atoms to consider>\n";
  usage += "\t@ref <reference coordinates>\n";
  usage += "\t@traj <trajectory files>\n";
  usage += "\t@step <read every n-th trajectory frame>\n";
  usage += "\t@ctraj <clustering trajectory file>\n";
  usage += "\t@normsd\n";
  usage += "\t@mat <clustering matrix file>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // get simulation time
    double time=0, dt=1;
    {
      Arguments::const_iterator iter=args.lower_bound("time");
      if(iter!=args.upper_bound("time")){
	time=atof(iter->second.c_str());
	++iter;
      }
      if(iter!=args.upper_bound("time"))
	dt=atof(iter->second.c_str());
    }

    // get normsd argument
    bool dormsd=true;
    try{
      args.check("normsd",0);
      dormsd=false;
    }
    catch(const Arguments::Exception &){
      dormsd=true;
    }

    // get step argument
    int frame_step=1;
    { 
      Arguments::const_iterator iter=args.lower_bound("step");
      if(iter!=args.upper_bound("step")){
        frame_step=atoi(iter->second.c_str());
      }
    }
    if (frame_step <= 0) {
      char message_buffer[120] = "";
      ostrstream ErMs(message_buffer,100);
      ErMs << "  Invalid value for step: " << frame_step << endl;
      throw gromos::Exception("cluster", message_buffer );
    }

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

    // open matrix file
    ofstream mf;
    string mfname_str="";
    char mfname[255];
    if (dormsd) {
      try{
        args.check("mat",1);
        mfname_str=""+args["mat"]+" ";
        mfname_str[mfname_str.length()-1]='\0';
        for (unsigned int i=0; i<mfname_str.length(); i++) 
          mfname[i]=mfname_str[i];
        mf.open(mfname);
      }
      catch(const Arguments::Exception &){
        mfname_str="./cluster_matrix_file.dat ";
        mfname_str[mfname_str.length()-1]='\0';
        for (unsigned int i=0; i<mfname_str.length(); i++)
          mfname[i]=mfname_str[i];
        mf.open(mfname);
      }
      if (!mf) {
        char message_buffer[120] = "";
        ostrstream ErMs(message_buffer,100);
        ErMs << "  Error opening matrix file " << mfname_str << endl;
        throw gromos::Exception("cluster", message_buffer );
      }
    }

    // open output trajectory file
    ofstream otf;
    OutG96 ot;
    string otfname_str="";
    char otfname[255];
    try{
      args.check("otraj",1);
      otfname_str=""+args["otraj"]+" ";
      otfname_str[otfname_str.length()-1]='\0';
      for (unsigned int i=0; i<otfname_str.length(); i++) 
        otfname[i]=otfname_str[i];
      otf.open(otfname);
      ot.open(otf);
      ot.writeTitle("Trajectory output file from clustering program");
    }
    catch(const Arguments::Exception &){
      otfname_str="./cluster_trajectory_file.dat ";
      otfname_str[otfname_str.length()-1]='\0';
      for (unsigned int i=0; i<otfname_str.length(); i++)
        otfname[i]=otfname_str[i];
      otf.open(otfname);
      ot.open(otf);
      ot.writeTitle("Trajectory output file from clustering program");
    }
    if (!otf) {
      char message_buffer[120] = "";
      ostrstream ErMs(message_buffer,100);
      ErMs << "  Error opening output trajectory file " << otfname_str << endl;
      throw gromos::Exception("cluster", message_buffer );
    }

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // gather reference system
    pbc->gather();
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
    // add classes
    for(Arguments::const_iterator it=args.lower_bound("class");
	it != args.upper_bound("class"); ++it){
      for(vector<int>::const_iterator mol=mols.begin();
	  mol!=mols.end();++mol)
	ref.addClass(*mol,it->second);
      added=1;
    }
    // add single atoms
    for(Arguments::const_iterator it=args.lower_bound("atoms");
	it != args.upper_bound("atoms"); ++it){
      int atom=atoi(it->second.c_str())-1, mol=0;
      while(atom >= refSys.mol(mol).numAtoms()){
	atom-=refSys.mol(mol).numAtoms();
	++mol;
	if(mol==refSys.numMolecules())
	  throw Arguments::Exception(usage);
      }
      ref.addAtom(mol,atom);
      added=1;
    }
    // did we add anything at all?
    if(!added)
      throw Arguments::Exception(usage);

    // System for calculation
    System sys(refSys);

    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);

    RotationalFit rf(&ref);
    Rmsd rmsd(&ref);

    // Store coordinates of first molecule in trajectories in array
    bool status=true;
    unsigned int frame_number=0;
    TrajArray ta(sys,MAXTRAJSIZ);

    // loop over all trajectories
    int frame_step_count=frame_step-1;

    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
       
      // loop over all frames
      while(!ic.eof()){
	ic >> sys; frame_step_count++;

        if (frame_step_count == frame_step) {
          frame_step_count=0;
	  pbc->gather();

	  rf.fit(&sys);

          // copy coordinates from sys and store
          status=ta.store(sys,frame_number);
          if (status) {
            frame_number++;
          }
          else {
            char message_buffer[120] = "";
            ostrstream ErMs(message_buffer,100);
            ErMs << "  Error storing frame number " << frame_number << endl;
            throw gromos::Exception("cluster", message_buffer );
          }

          // write to output trajectory
          ot << sys;
        }
      }
      ic.close();
    }
    ot.close();

    // Fit and calculate rmsd for all pairs with saved trajectory data 
    if (dormsd) {
      for (unsigned int stframe=0; stframe< frame_number-1 ; stframe++) {

        // extract stored coordinates and copy back into ref
        status=ta.extract(ref.sys(),stframe);
        if (!(status)) {
          char message_buffer[120] = "";
          ostrstream ErMs(message_buffer,100);
          ErMs << "  Error extracting reference frame " << stframe << endl;
          throw gromos::Exception("cluster", message_buffer );
        }

        for (unsigned int frame=stframe+1; frame < frame_number; frame++) {
          // extract stored coordinates and copy back into sys
          status=ta.extract(sys,frame);
          if (!(status)) {
            char message_buffer[120] = "";
            ostrstream ErMs(message_buffer,100);
            ErMs  << "  Error extracting frame " << frame << endl;
            throw gromos::Exception("cluster", message_buffer );
          }
  
  
          rf.fit(&sys);

          double r = rmsd.rmsd(sys);

          mf << setw(5) << stframe;
          mf << setw(5) << frame;
          mf.precision(5);
          mf << setw(11) << r;
          mf << endl;
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

