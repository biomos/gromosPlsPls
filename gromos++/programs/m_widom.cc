//ener calculates (non-bonded) interaction energies for specific atoms

#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/PtTopology.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InPtTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/utils/Energy.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/physics.h"
#include "../src/gmath/Distribution.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

  
int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "intopo", "insx", "inpttopo", 
		    "time", "stride", "cut", "eps", "kap",
                    "rdf", "rdfparam", "temp", "ntry", "traj"};
  int nknowns = 15;

  string usage = argv[0];
  usage += "\n\t@topo     <topology>\n";
  usage += "\t@pbc      <boundary type>\n";
  usage += "\t@intopo   <topology of the inserted particle>\n";
  usage += "\t@inpttopo <(multiple) perturbation topology>\n";
  usage += "\t@insx     <coordinates of the inserted particle>\n";
  usage += "\t@time     <time> <dt>\n";
  usage += "\t@stride   <take every n-th frame>\n";
  usage += "\t@cut      <cut-off distance>\n";
  usage += "\t@eps      <epsilon for RF>\n";
  usage += "\t@kap      <kappa for RF>\n";
  usage += "\t@rdf      <rdf with atom types>\n";
  usage += "\t@rdfparam <rdf-cutoff> <grid>\n";
  usage += "\t@temp     <temperature>\n";
  usage += "\t@ntry     <number of insertion tries per frame>\n";
  usage += "\t@traj     <trajectory files>\n";
  
 
try{
  clock();
  
  Arguments args(argc, argv, nknowns, knowns, usage);

  //   get simulation time
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

  // get the stride
  int stride=1;
  if(args.count("stride")>0) stride=atoi(args["stride"].c_str());

  //  read topologies, create a system that will contain only the topology.
  InTopology it(args["topo"]);
  System systop(it.system());
  GromosForceField gff(it.forceField());

  // count the number of atoms in the system
  int start=0;
  for(int m=0; m< systop.numMolecules(); ++m) start+=systop.mol(m).numAtoms();
  
  InTopology it2(args["intopo"]);
  System insys(it2.system());

  InPtTopology ipt(args["inpttopo"], start);
  PtTopology pert(ipt.ptTopo());

  // define input coordinate
  InG96 ic;
  ic.open(args["insx"]);
  ic >> insys;
  ic.close();
  
  // move the insertion particle to its first atom
  // usually it will be only one atom, but we keep our possibilities
  // open
  for(int m=0; m< insys.numMolecules(); m++)
    for(int a=0; a< insys.mol(m).numAtoms(); a++)
      insys.mol(m).pos(a)-=insys.mol(0).pos(0);
  
  // get non-bonded parameters
  double cut=atof(args["cut"].c_str());
  double eps=1.0, kap=0.0;
  if(args.count("eps")>0) eps=atof(args["eps"].c_str());
  if(args.count("kap")>0) kap=atof(args["kap"].c_str());
  
  // which atom types
  vector<AtomSpecifier> rdfatoms;
  {
    Arguments::const_iterator iter=args.lower_bound("rdf");
    Arguments::const_iterator to=args.upper_bound("rdf");
    for(int i=-1;iter!=to; ++iter){
      if(iter->second =="new"){
	AtomSpecifier as(systop);
	rdfatoms.push_back(as);
	i++;
      }
      else{
	rdfatoms[i].addSpecifier(iter->second);
      }
    }
  }
  
  // prepare rdf arrays
  double rdfcut=cut;
  int rdfgrid=100;
  {
    Arguments::const_iterator iter=args.lower_bound("rdfparam");
    if(iter!=args.upper_bound("rdfparam")){
      rdfcut=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("rdfparam"))
      rdfgrid=atoi(iter->second.c_str());
  }
  vector <gmath::Distribution> rdf;
  vector <vector<double> > s_rdf;

  for(unsigned int i=0; i<rdfatoms.size(); i++){
    gmath::Distribution d(0, rdfcut, rdfgrid);
    rdf.push_back(d);
  }
  s_rdf.resize(rdfatoms.size()*pert.numPt());
  for(unsigned int i=0; i< s_rdf.size(); i++){
    s_rdf[i].resize(rdfgrid, 0.0);
  }
  
  double rdfdist=rdfcut/double(rdfgrid);
  double rdfcorr=1/(4.0*acos(-1.0)*rdfdist);
  
  // read in the temperature
  double temp=atof(args["temp"].c_str());
  double beta=-1.0/BOLTZ/temp;
  
  // read in number of tries per frame
  int ntry=atoi(args["ntry"].c_str());
  
  // define some values for averaging
  int numframes=0;
  double fexp=0;               // exp(-E/kT)
  
  double vol, s_vol=0.0;       // volume
  double v_exp;
  vector<double> s_v_exp(pert.numPt(), 0.0);   // v.exp(-E/kt)
  vector<double> s_v_Eexp(pert.numPt(),0.0);   // v.E.exp(-E/kt)
  
  // initialize random seed
  srand(int(clock()));
  cout << setw(6) << "# time"
       << setw(22) << "DG"
       << setw(22) << "DH(solu-solv)"
       << setw(22) << "TDS (DH-DG)"
       << endl;

  // loop over all trajectories
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
    ic.open((iter->second).c_str());
    ic.select("ALL");
    
      // loop over single trajectory
    while(!ic.eof()){
      numframes++;

      // create a new system to which we can add the molecules of the 
      // insys;
      System sys(systop);
      
      // parse boundary conditions

      Boundary *pbc = BoundaryParser::boundary(sys, args);
      Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
      
      // read in the coordinates
      ic >> sys;
      
      if(!(numframes % stride)){
	// we have to gather to get covalent interactions 
	// and charge-groups connected

	(*pbc.*gathmethod)();
	
	// add the molecules of the insys
	for(int m=0; m<insys.numMolecules(); ++m){
	  sys.addMolecule(insys.mol(m));
	  sys.mol(systop.numMolecules() + m).initPos();
	}
	
	// create the energy class for the new system
	Energy en(sys, gff, *pbc);
	
	// set the atoms for which to calculate the energy
	AtomSpecifier energyatoms(sys);
	for(int m=0; m<insys.numMolecules(); ++m){
	  for(int a=0; a<insys.mol(m).numAtoms(); ++a){
	    energyatoms.addAtom(systop.numMolecules() + m, a);
	  }
	}
	en.setAtoms(energyatoms);
	en.setCutOff(cut);
	en.setRF(eps, kap);
	
	// now reset the rdfatoms to point at this system
	for(unsigned int i=0; i< rdfatoms.size(); ++i){
	  rdfatoms[i].setSystem(sys);
	}
	
	// we need the volume, correct for truncated octahedron!!
	vol=sys.box()[0]*sys.box()[1]*sys.box()[2];
	if(pbc->type()=='t') vol /= 2;
	s_vol+=vol;
	
	// now we can go into the loop over the trial positions
	for(int i=0; i<ntry; i++){
	  //get three random numbers between the box dimensions
	  Vec move;
	  
	  for(int d=0; d<3; d++){
	    int r=rand();
	    move[d]=double(sys.box()[d]*r)/double(RAND_MAX);
	  }
	  
	  // move the inatoms
	  for(int m=0; m< insys.numMolecules(); ++m){
	    for(int a=0; a< insys.mol(m).numAtoms(); ++a){
	      sys.mol(systop.numMolecules() + m ).pos(a) = 
		insys.mol(m).pos(a) + move;
	    }
	  }
	  
	  // do the rdf's
	  for(unsigned int j=0; j<rdf.size(); j++){
	    
	    // set them to zero
	    rdf[j].clear();
	    
	    // loop over the 'with' atoms
	    for(int i=0; i < rdfatoms[j].size(); i++){
	      Vec tmp;
	      tmp=pbc->nearestImage(move, *rdfatoms[j].coord(i),
				    sys.box());
	      rdf[j].add((tmp-move).abs());
	    }
	  }

	  // For the actual energy calculation create the pairlist
	  en.calcPairlist();

	  // loop over the different perturbations
	  for(int p=0; p < pert.numPt(); ++p){
	    // set the parameters
	    pert.apply(sys, p);
	    
	    // calculate the interactions
	    en.calcNb_interactions();
	    
	    // store and sum everything to the appropriate arrays
	    fexp=exp(beta*en.tot());
	    v_exp=vol*fexp;
	    
	    s_v_exp[p]+=v_exp;
	    s_v_Eexp[p]+=en.tot()*v_exp;
	    
	    for(unsigned int j=0; j<rdfatoms.size(); j++){
	      for(int i=0; i< rdfgrid; i++){
		double r=(i+0.5)*rdfdist;
		
		s_rdf[p*rdfatoms.size() + j][i] += 
		  v_exp*rdf[j][i]*vol*rdfcorr/(r*r*rdfatoms[j].size());
	      }
	    }
	  } // loop over perturbations
	} // loop over trials
	cout << setw(6) << time;
	for(int p=0; p<pert.numPt(); p++){
	  double DG=-BOLTZ*temp*log(s_v_exp[p]/s_vol/ntry);
	  double DH= s_v_Eexp[p]/s_v_exp[p];
	  double TDS=(DH-DG);
	  
	  cout << setw(12) << DG
	       << setw(12) << DH
	       << setw(12) << TDS;
	}
	cout << endl;
      } // if stride



      
      time+=dt;
    } //loop over frames
  }
  // print out averages

  for(int p=0; p< pert.numPt(); ++p){
    cout << "# ---------------------\n"
	 << "# " << pert.pertName(p) << endl;
  
    cout.precision(10);
    cout.setf(ios::right, ios::adjustfield);
    double DG=-BOLTZ*temp*log(s_v_exp[p]/s_vol/ntry);
    double DH=s_v_Eexp[p]/s_v_exp[p];
    double DS=(DH-DG)/temp;
    
    cout << endl;
    cout << "# Number of frames: " << numframes << endl;
    cout << "# Number of tries per frame: " << ntry << endl;
    cout << "# <V.exp(-E/kT)> : " << s_v_exp[p]/numframes/ntry << endl;
    cout << "# <V.E.exp(-E/kT)> : " << s_v_Eexp[p]/numframes/ntry << endl;
    cout << "# <V> : " << s_vol/numframes << endl;
    cout << "# DG: " << DG << endl;
    cout << "# DH: " << DH << endl;
    cout << "# DS: " << DS << endl;
    

    ostringstream os;
    os << "rdf_widom_" << p+1 << ".out";
  
    ofstream fout(os.str().c_str());
    fout << "# rdf's for insertion of:" << endl;
    fout << "#    " << args["insx"] << endl;
    fout << "# into:" << endl;
    fout << "#    " << args["traj"] << endl;
    fout << "#" << endl;
    fout << "# taking rdf of test position with "
	 << rdfatoms.size() << " specified groups of atoms\n";
    fout << "#" << endl;
    fout << "#" << setw(11) << "r";
    for(unsigned int j=0; j< rdfatoms.size(); j++)
      fout << " " << setw(8) << j;
    fout << endl;
    for(int i=0; i<rdfgrid; i++){
      
      fout << setw(12) << (i+0.5)*rdfdist;
      for(unsigned int j=0; j< rdfatoms.size(); j++)
	fout << " " << setw(12) << s_rdf[p*rdfatoms.size() + j][i]/s_v_exp[p];
      fout << endl;
      
    }
    fout.close();
  } // loop over perturbations
  
}
 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}









