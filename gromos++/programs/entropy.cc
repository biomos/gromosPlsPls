// Version with ROT and TRANS contributions REMOVED.
// -> There is always fitting and putting solute into center of box
// remove all solvent options
// can do Schlitter and Quasiharmonic analysis

#include <cassert>
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/Rmsd.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/PositionUtils.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>

using namespace utils;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace fit;
using namespace std;

// Function declaration ---------------------------------------------

// computing the entropy (Schlitter approximation)
double entropy(gsl_matrix * cov, gsl_vector * av, gsl_vector * mass, unsigned long confs, int ndof, double temperature);
// computing the entropy with the quasiharmonic analysis
double freq_etc(gsl_matrix * cov,gsl_vector * av,gsl_vector * mass,unsigned long confs,int ndof,double temperature,int printevs);

int readRef(System &refsys, Reference &ref, Arguments &args, string usage);
void readEatoms(System & sys, AtomSpecifier & entropy_atoms, Arguments & args, string usage);

void print_setting(int nat_solu, AtomSpecifier & at_S, System &sys, gsl_vector *mass, int fit, Reference &ref, string how = "nv");

// define necessary constants
const double KB=1.38065e-23;        // Boltzmans constant
const double E=2.718281828459;      // Euler number
const double HBAR=1.054571596e-34;  // Planks constant over 2 Pi
const double MU=1.66053873e-27;     // Atomic mass unit
const double NA=6.022142e23;        // Avogadros number

int main(int argc, char *argv[]){

  char *knowns[] = {"topo", "pbc", "pbc_ref",
		    "fit", "atoms_fit", "atoms_entropy", 
		    "temp",
		    "win_ave", "time", 
		    "meth" , "traj"};

  int nknowns = 12;
  
  string usage = argv[0];
  usage += "\n\t@topo           <topology>\n";
  usage += "\t@pbc            <boundary type>\n";
  usage += "\t@fit            <structure to fit against>\n";
  usage += "\t# @pbc_ref      <boundary type for reference for fit>\n";
  usage += "\t@atoms_fit      <atoms to consider for fit (atomspecifier)>\n";
  usage += "\t@atoms_entropy  <atoms to consider for Entropy (atomspecifier)>\n";
  usage += "\t@temp           <temperature>\n";
  usage += "\t@time           <dt stride blocklength [skip stop]>\n";
  usage += "\t# @win_ave      <averaging over windows of length nana (1 or 0 / default:0)>\n";
  usage += "\t# @meth         <methods to use schlitter=yes, quasiharm=yes (default both)>\n";
  usage += "\t@traj           <trajectory files>\n";
   
  try{
      Arguments args(argc, argv, nknowns, knowns, usage);

      //  read topology -------------------------------
      args.check("topo",1);
      InTopology it(args["topo"]);
      System sys(it.system());

      // parse the type of boundary conditions and create pbc
      Boundary *pbc = BoundaryParser::boundary(sys, args, "pbc"); 
     // parse the gathering method
      Boundary::MemPtr gathmethod = args::GatherParser::parse(args, "pbc");
          
      //   get simulation temperature ------------------
      double temp; 
      temp=atof(args["temp"].c_str());

      //   get simulation time-step dt [ps]
      //   the trajectory is read every nread steps
      //   intervals for analysis: nana steps (0 -> only at the end)
      double dt=1; 
      unsigned long nread=1,nana=0,start=0,end=0;
      {
	  Arguments::const_iterator iter=args.lower_bound("time");
	  if(iter!=args.upper_bound("time")){
	    dt=atof(iter->second.c_str());
	    ++iter;
	  }
	  if(iter!=args.upper_bound("time")){
	    nread=atoi(iter->second.c_str());
	    ++iter;
	  }  
	  if(iter!=args.upper_bound("time")){
	    nana=atoi(iter->second.c_str());
	    ++iter;
	  }
	  if(iter!=args.upper_bound("time")){
	    start=atoi(iter->second.c_str());
	    ++iter;
	  }
	  if(iter!=args.upper_bound("time"))
	    end=atoi(iter->second.c_str());
      }

      //   averaging over windows / yes or no ?  ------------------
      int window_averaging = 0; 
      {
	  Arguments::const_iterator iter=args.lower_bound("win_ave");
	  if(iter!=args.upper_bound("win_ave")){
	      window_averaging=atoi(iter->second.c_str());
	  }
      }
   
      //  which method: schlitter, quasiharmonic, both -------------
      int schlitter=1, quasi=1; 
      {
	  Arguments::const_iterator iter=args.lower_bound("meth");
	  if(iter!=args.upper_bound("meth")){
	    schlitter=atoi(iter->second.c_str());
	    ++iter;
	  }
	  if(iter!=args.upper_bound("meth")){
	    quasi=atoi(iter->second.c_str());
	  }
      }
      if (!(schlitter || quasi)) {	
	cerr << "You should do either schlitter or quasi!!!!\n";
	return -1;
      }

      // only solute molecules -------------------------
      AtomSpecifier atoms_entropy(sys);
      readEatoms(sys, atoms_entropy, args, usage);
      
      const int nat_solu = atoms_entropy.size();
      const int ndof = 3 * nat_solu; // nr of degrees of freedom for entropy

      // fitting -------------------------------------------
      System refsys(sys);
      Reference ref(&refsys);      
      int fit = readRef(refsys,ref,args,usage);
      
      // define the necessary vectors and matrices ----------------
      // averages 
      gsl_vector *pos_av = gsl_vector_alloc (ndof);
      gsl_vector_set_zero (pos_av);
      // position vector to store coordinates intermediately
      gsl_vector *position = gsl_vector_alloc (ndof);
      // masses (per degree of freedom, easier to calculate...)
      gsl_vector *mass = gsl_vector_alloc(ndof);

      // matrices
      gsl_matrix * covariance = gsl_matrix_alloc (ndof, ndof);
      gsl_matrix_set_zero (covariance);
      
      // fill the mass vector
      // easy to replace with atom spec...
      int count = 0;
      for(int i=0; i<atoms_entropy.size(); ++i){
	for(int k=0; k<3; ++k){

	  gsl_vector_set(mass, count, atoms_entropy.mass(i));
	  ++count;
	}
      }
      
      
      // print the setting of all parameters (adapt!)
      print_setting(nat_solu, atoms_entropy, sys, mass, fit, ref, "nv");
      // temperature etc ...
    
      // define input coordinate
      InG96 ic;
      RotationalFit rf(&ref);
      unsigned long numFrames = 0;
      unsigned long allFrames = 0;
      

      // loop over all trajectories
      for(Arguments::const_iterator iter=args.lower_bound("traj"),
	    to=args.upper_bound("traj");iter!=to; ++iter){

	ic.open((iter->second).c_str());
	ic.select("SOLUTE");
      
	while(!ic.eof()){

	  ic >> sys;	  
	  ++allFrames;
		
	  if((allFrames > start) && ((allFrames - start - 1) % nread == 0)) {   
	    ++numFrames;

	    (*pbc.*gathmethod)();

	    if (fit) {
	      try {
		rf.fit(&sys);
	      }
	      catch (RotationalFit::Exception &degenerate){
		cerr << allFrames;
		cerr << ": this Frame is analysed without fitting (maybe try Kabsch-fitting)\n";
	      }
	    }

	    // get the positions
	    // replace by atom spec...
	    count = 0;
	    for(int i=0; i<atoms_entropy.size(); ++i){
	      for(int k=0; k<3; ++k){

		gsl_vector_set(position, count, atoms_entropy.pos(i)[k]);
		count++;
	      }
	    }
	  
	    // fill the average vector and covariance matrix
	    for(int i=0; i< ndof; ++i){
	      double pos_i;
	      pos_i = gsl_vector_get(position, i);
	      gsl_vector_set(pos_av, i, gsl_vector_get(pos_av, i) + pos_i);
	      for(int j=i; j < ndof; ++j)
		gsl_matrix_set(covariance, i, j,
			       (gsl_matrix_get(covariance, i, j) + pos_i * gsl_vector_get(position, j))
			       );
	    }
	    
	    // (cumulative) entropy time-series
	    if ((nana > 0) && (numFrames > 1) &&
		((allFrames-start-1) % nana == 0)){

	      cout << (allFrames-1)*dt << "  ";
	      if (schlitter)
		cout << entropy(covariance, pos_av, mass, numFrames,
				ndof, temp);
	      cout << "   ";
	      if (quasi)
		cout  << freq_etc(covariance, pos_av, mass, numFrames,
				  ndof, temp, 0);
	      cout << endl;
	      // reset to zero
	      if (window_averaging) {
		gsl_vector_set_zero (pos_av);
		gsl_matrix_set_zero (covariance);
		numFrames = 0;
	      }
	    }	    
	  }
	  if (allFrames == end) break;

	} // loop over frames in trajectory
	ic.close();   

	if (allFrames == end) break;
      } // loop over trajectories

      if (! window_averaging &&
	  !((allFrames - start - 1) % nana == 0)){

	cout  << (allFrames-1) * dt << "  ";;
	if (schlitter)
	  cout << entropy(covariance, pos_av, mass, numFrames, ndof, temp) ;
	cout << "  ";
	if (quasi)
	  cout  << freq_etc(covariance,pos_av,mass,numFrames,ndof,temp,1) ;
	cout << endl;
      }

      cout << "\n# entropy was analysed using " << numFrames - start << " frames\n";

      gsl_matrix_free(covariance);   
      gsl_vector_free(mass);
      gsl_vector_free(pos_av);
      gsl_vector_free(position);
      
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    return 1;
  }
  return 0;
}


// Function to compute the entropy ----------------------------------
double entropy(gsl_matrix * cov, gsl_vector * av, gsl_vector * mass, 
	       unsigned long confs, int ndof, double temperature) {
  
  // confs : number of frames in average positions and covariance matrix

  // calculate factor: mu*k*T*e^2/hbar^2 (+conversion from nm to m)
  double mukte2h2 = temperature * 1e-18 * MU * KB * E * E / (HBAR * HBAR);

  gsl_matrix * lu = gsl_matrix_alloc (ndof, ndof);

  // diagonal elements
  for (int i=0; i<ndof; ++i) {
    double lui = gsl_vector_get(av, i) / confs;
    lui = - lui * lui + gsl_matrix_get(cov, i, i) / confs;
    lui = 1 + lui * mukte2h2 * gsl_vector_get(mass, i);
    gsl_matrix_set(lu, i, i, lui);
  }
  // off-diagonal elements
  for (int i=0; i<ndof; ++i) {
    for (int j=i+1; j<ndof; ++j) {
      double luij = gsl_matrix_get(cov, i, j) / confs - gsl_vector_get(av, i) * gsl_vector_get(av, j) / (confs * confs);
      luij *= mukte2h2;
      gsl_matrix_set(lu, i, j, luij * gsl_vector_get(mass,j));
      gsl_matrix_set(lu, j, i, luij * gsl_vector_get(mass,i));
    }
  }

  gsl_permutation * p = gsl_permutation_alloc(ndof);
  int s;
  gsl_linalg_LU_decomp(lu, p, &s);

  double lndet = gsl_linalg_LU_lndet(lu);
  double entropy = 0.5 * KB * NA * lndet;
   
  gsl_permutation_free(p);
  gsl_matrix_free(lu);

  return entropy;
}


double freq_etc(gsl_matrix * cov, gsl_vector * av, gsl_vector * mass, 
		unsigned long confs, int ndof, double temperature, int printevs) {

  // calculate factor: mu*k*T*e^2/hbar^2 (+conversion from nm to m)
  // double mukte2h2 = temperature * 1e-18 * MU * KB * E * E / (HBAR * HBAR);
  double mufac =  1e-18 * MU;
  double kT = temperature * KB;
  
  gsl_matrix * lu = gsl_matrix_alloc (ndof, ndof);

  // diagonal elements
  for (int i=0; i < ndof; ++i) {
    double lui = gsl_vector_get(av, i) / confs;
    lui = - lui * lui + gsl_matrix_get(cov, i, i) / confs;
    lui *= mufac * gsl_vector_get(mass, i);
    gsl_matrix_set(lu, i, i, lui);
  }
  // off-diagonal elements
  for (int i=0; i < ndof; ++i) {
    for (int j=i+1; j < ndof; ++j) {
      double luij = gsl_matrix_get(cov, i, j) / confs -
	gsl_vector_get(av, i) * gsl_vector_get(av, j) / (confs * confs);
      luij *= mufac;
      // for a faster program -> use a sqrtmass_vector
      luij *= sqrt(gsl_vector_get(mass, i)) * sqrt(gsl_vector_get(mass, j));      
      // MARKUS: this one is now symmetric, not like the above one??
      gsl_matrix_set(lu, i, j, luij);
      gsl_matrix_set(lu, j, i, luij);
    }
  }

  gsl_vector *eval = gsl_vector_alloc (ndof);
  gsl_matrix *evec = gsl_matrix_alloc (ndof, ndof);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (ndof);
  gsl_eigen_symmv (lu, eval, evec, w);
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  // MARKUS: this overwrites it for every output
  // (plus; the format seems a bit xmgrace - unfriendly)
  ofstream evfile; 
  if (printevs)
    evfile.open("Evals.out");

  double entropy = 0.0;
  double ho_kT;
  int nr_evals = 0;
  for (int i=0; i<ndof; ++i) {  
    double ev_i = sqrt(kT * gsl_vector_get(eval, i)) / HBAR;
    if (printevs)
      evfile << "# EV :  " << gsl_vector_get(eval, i) << "  ";
    if (ev_i > 1.0e-10){
      ++nr_evals;
      ho_kT = 1.0 / ev_i;
      ev_i = KB * NA * (ho_kT / (exp(ho_kT) - 1.0) - log(1.0 - exp(-ho_kT)));
      entropy += ev_i;
      if (printevs)
	evfile << "  -  contrib to entropy: " << ev_i;
    }
    if (printevs)
      evfile << endl;    
  }
  if (printevs)  
    evfile << "# EV : " << nr_evals << " nr of eigenvalues.\n";

  if (printevs)
    evfile.close();
  
  gsl_eigen_symmv_free(w); 
  gsl_matrix_free(lu);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);

  return entropy;
}

int readRef(System &refsys, Reference &ref, Arguments &args, string usage) 
{
  if (args.count("fit") >= 0){

    InG96 fitRef;
    fitRef.open(args["fit"]);
    fitRef.select("SOLUTE");
    fitRef >> refsys;

    // PBC
    if (args.count("pbc_ref") > 0) {
      Boundary *refpbc = BoundaryParser::boundary(refsys, args, "pbc_ref"); 
      Boundary::MemPtr refgathmethod = args::GatherParser::parse(args, "pbc_ref");
      (*refpbc.*refgathmethod)();
      delete refpbc;
    }
    else{ // take standard "pbc"
      Boundary *refpbc = BoundaryParser::boundary(refsys, args); 
      Boundary::MemPtr refgathmethod = args::GatherParser::parse(args);
      (*refpbc.*refgathmethod)();
      delete refpbc;
    }
    
    fitRef.close();

    // get atoms for fit
    AtomSpecifier atoms_fit(refsys);
    Arguments::const_iterator iter = args.lower_bound("atoms_fit");
    Arguments::const_iterator to = args.upper_bound("atoms_fit");
    
    for(; iter != to; iter++){
      atoms_fit.addSpecifier(iter->second);
    }

    ref.addAtomSpecifier(atoms_fit);
    ref.normalise();
    
    return 1;
  }
	 
  return 0;
}


void print_setting(int nat_solu, AtomSpecifier & at_S, System &sys,
		   gsl_vector *mass, int fit, Reference &ref, 
		   string how)
{
       
  cout << "# NAT solute: " << nat_solu << endl;
  
  int count = 0;
  for(int i=0; i<at_S.size(); ++i){
    if (how == "v") {  
      cout << "# mol: " << at_S.mol(i) << " " << "at: " << at_S.atom(i) << " ";
      cout << at_S.name(i);
      cout << "  mass: ";
      cout << at_S.mass(i) << endl;
    }
  }

  if (fit){
    count = 0;
    if (how == "v") 
      cout << "# Fitting is done on the following atoms: \n";
    for (int i=0;i<sys.numMolecules();++i)
      for(int j=0;j < sys.mol(i).numAtoms();++j)
	if (ref.weight(i,j) > 0) {
	  if (how == "v") {
	    cout << "# mol: " << i << " " << "at: " << j << " ";
	    cout << sys.mol(i).topology().atom(j).name();
	    cout << "  weight: " << ref.weight(i,j) << endl;
	  }
	  count += 1;
	}
    cout << "# Fitting is done on " << count << " atoms\n";
  } else cout << "# No fitting\n";
  
  return;
}

void readEatoms(System & sys, AtomSpecifier & entropy_atoms, Arguments & args, string usage) 
{
  Arguments::const_iterator iter=args.lower_bound("atoms_entropy");
  Arguments::const_iterator to=args.upper_bound("atoms_entropy");

  for( ; iter!=to; ++iter){
    entropy_atoms.addSpecifier(iter->second);
  } 
}
