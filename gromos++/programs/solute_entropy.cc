/**
 * @file solute_entropy.cc
 * Calculates Schlitter and Quasiharmonic analysis
 */
/**
 * @page programs Program Documentation
 *
 * @anchor solute_entropy
 * @section solute_entropy Schlitter and Quasiharmonic analysis
 * @author @ref rba @ref ns
 * @date 15-09-08
 *
 * Program solute_entropy takes a coordinate trajectory and calculated the
 * configurational entropy using the Schlitter and the quasiharmonic analysis 
 * methods for a given set of atoms.
 * The entropy can be averaged over a window of a given size.
 * If requested, a rotational fit prior to the entropy calculation is carried out.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;time and dt&gt; </td></tr>
 * <tr><td> \@atomsentropy</td><td>&lt;atomspecifier: atoms to consider for entropy&gt; </td></tr>
 * <tr><td> \@temp</td><td>&lt;temperature&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td>[\@ref</td><td>&lt;reference structure to fit against&gt;]</td></tr>
 * <tr><td>[\@ref_pbc</td><td>&lt;boundary type for reference for fit&gt;]</td></tr>
 * <tr><td>[\@atomsfit</td><td>&lt;atomspecifier: atoms to consider for fit&gt;]</td></tr>
 * <tr><td>[\@average</td><td>&lt;averaging over windows of this length&gt;]</td></tr>
 * <tr><td>[\@method</td><td>&lt;methods to use: schlitter quasiharm (default both)&gt;]</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  solute_entropy
    @topo       ex.top
    @pbc        r
    @time       0 0.1
    @atomsentropy  1:CB
    @atomsfit   1:CB
    @ref        ex.coo
    @ref_pbc    r
    @traj       ex.tr
    @temp       300
    @method     schlitter quasiarm
    @average    100
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gio/OutPdb.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/physics.h"
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
// Function to compute the entropy ----------------------------------
double entropy(gsl_matrix * cov, gsl_vector * av, gsl_vector * mass, unsigned long confs, int ndof, double temperature);
// computing the entropy with the quasiharmonic analysis
double freq_etc(gsl_matrix * cov, gsl_vector * av, gsl_vector * mass, unsigned long confs, int ndof, double temperature);


// define necessary constants
const double KB = 1.38065e-23; /// Boltzmans constant
const double E = 2.718281828459; /// Euler number
const double HBAR = 1.054571596e-34; /// Planks constant over 2 Pi
const double MU = AMU; /// Atomic mass unit
const double NA = AVOGADRO; /// Avogadros number

int main(int argc, char **argv) {

  char *knowns[] = {"topo", "pbc", "ref" ,"ref_pbc", "atomsfit", "temp", "atomsentropy", "average", "time", "method", "traj"};
  int nknowns = 11;

  string usage = "# ";
  usage += string(argv[0]) + "\n";
  usage += "\n\t@topo        <topology>\n";
  usage += "\t@pbc           <boundary type>\n";
  usage += "\t@ref           <structure to fit against>\n";
  usage += "\t@ref_pbc       <boundary type for reference for fit>\n"; 
  usage += "\t@atomsfit      <atoms to consider for fit (atomspecifier)>\n";
  usage += "\t@atomsentropy  <atoms to consider for Entropy (atomspecifier)>\n"; 
  usage += "\t@temp          <temperature>\n";
  usage += "\t@time          <time and dt>\n";
  usage += "\t@average       <averaging over windows of this length>\n";
  usage += "\t@method        <methods to use schlitter=yes, quasiharm=yes (default both)>\n";
  usage += "\t@traj          <trajectory files>\n";

  try {
    Arguments args(argc, argv, nknowns, knowns, usage);


    //  read topology -------------------------------
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());

    // parse the type of boundary conditions and create pbc
    Boundary *pbc = BoundaryParser::boundary(sys, args, "pbc");
    // parse the gathering method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args, "pbc");

    //   get simulation temperature ------------------
    double temp = 0.0;
    {
      istringstream str(args["temp"]);
      if (!(str >> temp) || temp == 0.0)
        throw gromos::Exception("solute_entropy", "temperature is not numeric or zero.");
    }

    //   get simulation time-step dt [ps]
    //   the trajectory is read every nread steps
    //   intervals for analysis: nana steps (0 -> only at the end)
    double time = 0, dt = 1;
    {
      Arguments::const_iterator iter = args.lower_bound("time");
      if (iter != args.upper_bound("time")) {
        istringstream str(iter->second);
        ++iter;
        
        if (!(str >> time))
          throw gromos::Exception("solute_entropy", "@time time is not numeric.");
      }
      if (iter != args.upper_bound("time")) {
        istringstream str(iter->second);
        ++iter;
        
        if (!(str >> dt))
          throw gromos::Exception("solute_entropy", "@time dt is not numeric.");
      }
      if (iter != args.upper_bound("time")) {
        throw gromos::Exception("solute_entropy", "@time too many arguments");
      }
    }

    //   averaging over windows / yes or no ?  ------------------
    int window_averaging = 0;
    {
      Arguments::const_iterator iter = args.lower_bound("average");
      if (iter != args.upper_bound("average")) {
        istringstream str(iter->second);
        if (!(str >> window_averaging))
          throw gromos::Exception("solute_entropy", "@average is not an integer.");
      }
    }

    //  which method: schlitter, quasiharmonic, both -------------
    bool schlitter = true, quasi = true;
    {
      Arguments::const_iterator iter = args.lower_bound("method"),
      to = args.upper_bound("method");
      if (iter != to) {
        schlitter = false;
        quasi = false;
      }
      // search arguments for methods
      for(; iter != to; ++iter) {
        if (iter->second == "schlitter")
          schlitter = true;
        if (iter->second == "quasiharm")
          quasi = true;
        if (iter->second == "both")
          schlitter = quasi = true;
      }
    }
    if (!schlitter && !quasi) {
      throw gromos::Exception("solute_entropy", "@method You should do either schlitter or quasi.");
    }

    // only solute molecules -------------------------
    // use the reference class as container for the selected atoms
    Reference atoms_entropy(&sys);
    AtomSpecifier atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atomsentropy");
      Arguments::const_iterator to = args.upper_bound("atomsentropy");
      for (; iter != to; iter++) 
        atoms.addSpecifier(iter->second);
    }

    // add atoms to the reference
    atoms_entropy.addAtomSpecifier(atoms);
    if (atoms.empty()) 
      throw gromos::Exception("solute_entropy", "@atomsentropy no atoms specified.");

    int ndof = 3 * atoms.size(); // nr of degrees of freedom for entropy


    // fitting -------------------------------------------
    System refsys(sys);
    Reference ref(&refsys);
    bool fit = false;
    {
      Arguments::const_iterator iter = args.lower_bound("ref");
      if (iter != args.upper_bound("ref")) {
        fit = true;
        InG96 fitRef;
        fitRef.open(iter->second.c_str());
        fitRef.select("SOLUTE");
        fitRef >> refsys;
        fitRef.close();
        // correct for periodic boundary conditions by calling
        // the appropriate gathering method	      
        // parse the type of boundary conditions and create pbc
        if (args.count("ref_pbc") > 0) {
          Boundary *refpbc = BoundaryParser::boundary(refsys, args, "ref_pbc");
          // parse the gathering method
          Boundary::MemPtr refgathmethod = args::GatherParser::parse(args, "ref_pbc");
          (*refpbc.*refgathmethod)();
          delete refpbc;
        }
      }
    }

    if (fit) {
      AtomSpecifier ref_atoms(refsys);
      Arguments::const_iterator iter = args.lower_bound("atomsfit");
      Arguments::const_iterator to = args.upper_bound("atomsfit");
      for (; iter != to; iter++)
        ref_atoms.addSpecifier(iter->second);

      if (ref_atoms.empty())
        throw gromos::Exception("solute_entropy", "@atomsfit no atoms specified.");

      // add atoms to the reference
      ref.addAtomSpecifier(ref_atoms);
    }

    // define the necessary vectors and matrices ----------------
    // averages 
    gsl_vector *pos_av = gsl_vector_alloc(ndof);
    gsl_vector_set_zero(pos_av);
    // position vector to store coordinates intermediately
    gsl_vector *position = gsl_vector_alloc(ndof);
    // masses (seems to be 3 times too big 
    // -> leaves space for mass scaling of specific DOFs)
    gsl_vector *mass = gsl_vector_alloc(ndof);
    // matrices
    gsl_matrix * covariance = gsl_matrix_alloc(ndof, ndof);
    gsl_matrix_set_zero(covariance);


    // fill the mass vector
    int count = 0;
    for (int i = 0; i < sys.numMolecules(); ++i) {
      for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
        if (atoms_entropy.weight(i, j) > 0) {
          for (int k = 0; k < 3; ++k) {
            gsl_vector_set(mass, count, sys.mol(i).topology().atom(j).mass());
            count++;
          }
        }
      }
    }  
    
    cerr << "WARNING: this program is experimental." << endl;
    
    // define input coordinate
    InG96 ic;
    RotationalFit rf(&ref);
    unsigned long numFrames = 0;
    
    // loop over all trajectories
    for (Arguments::const_iterator iter = args.lower_bound("traj"),
            to = args.upper_bound("traj"); iter != to; ++iter) {

      ic.open((iter->second).c_str());

      while (!ic.eof()) {
        ic.select("SOLUTE");
        ic >> sys;

        // correct for periodic boundary conditions by calling
        // the appropriate gathering method	 
        (*pbc.*gathmethod)();
        if (fit)
          rf.fit(&sys);

        // print title
        if (numFrames == 0) {
          cout.setf(ios::right, ios::adjustfield);
          cout << "# " << setw(8) << "time" << setw(15) << "schlitter"
               << setw(15) << "quasi" << endl;
        }
        
        // get the positions
        count = 0;
        for (int i = 0; i < sys.numMolecules(); ++i) {
          for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
            if (atoms_entropy.weight(i, j) > 0) {
              for (int k = 0; k < 3; ++k) {
                gsl_vector_set(position, count, sys.mol(i).pos(j)[k]);
                count++;
              } // k
            } // if selected
          } // atoms
        } // molecules


        // fill the average vector/matrix
        double temporary;
        for (int i = 0; i < ndof; i++) {
          temporary = gsl_vector_get(position, i);
          gsl_vector_set(pos_av, i, gsl_vector_get(pos_av, i) + temporary);
          for (int j = i; j < ndof; j++)
            gsl_matrix_set(covariance, i, j, (gsl_matrix_get(covariance, i, j) + temporary * gsl_vector_get(position, j)));
        }

        cout.precision(2);
        cout.setf(ios::right, ios::adjustfield);
        cout << setw(10) << time;

        cout.precision(8);
        cout.setf(ios::fixed, ios::floatfield);
        
        double entr_schlitter = 0.0, entr_quasi = 0.0;
        
        if (schlitter)
          entr_schlitter = entropy(covariance, pos_av, mass, numFrames + 1, ndof, temp);
        
        if (quasi)
          entr_quasi = freq_etc(covariance, pos_av, mass, numFrames + 1, ndof, temp);
        
        cout << setw(15) << entr_schlitter
             << setw(15) << entr_quasi << endl;

        if (window_averaging && numFrames && numFrames % window_averaging == 0) {
          gsl_vector_set_zero(pos_av);
          gsl_matrix_set_zero(covariance);
        } // averaging
          
        time += dt;
        numFrames++;
      } // frames
      ic.close();
    } // files

    gsl_matrix_free(covariance);
    gsl_vector_free(mass);
    gsl_vector_free(pos_av);
    gsl_vector_free(position);

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

// Function to compute the entropy ----------------------------------
double entropy(gsl_matrix * cov, gsl_vector * av, gsl_vector * mass, unsigned long confs, int ndof, double temperature) {
  const double dbl_confs = double(confs);
  // calculate factor: mu*k*T*e^2/hbar^2 (+conversion from nm to m)
  double mukte2h2 = temperature * 1e-18 * MU * KB * E * E / (HBAR * HBAR);
  gsl_matrix * lu = gsl_matrix_alloc(ndof, ndof);

  double temp;
  // off-diagonal elements
  for (int i = 0; i < ndof; i++) {
    for (int j = i + 1; j < ndof; j++) {

      temp = gsl_matrix_get(cov, i, j) / dbl_confs - gsl_vector_get(av, i) * gsl_vector_get(av, j) / (dbl_confs * dbl_confs);
      temp = temp * mukte2h2;
      gsl_matrix_set(lu, i, j, temp * gsl_vector_get(mass, j));
      gsl_matrix_set(lu, j, i, temp * gsl_vector_get(mass, i));
    }
  }
  // diagonal elements
  for (int i = 0; i < ndof; i++) {
    temp = gsl_vector_get(av, i) / dbl_confs;
    temp = -temp * temp + gsl_matrix_get(cov, i, i) / dbl_confs;
    temp = 1 + temp * mukte2h2 * gsl_vector_get(mass, i);
    gsl_matrix_set(lu, i, i, temp);
  }


  double diag = 0.0;
  for (int i = 0; i < ndof; i++)
    diag += log(gsl_matrix_get(lu, i, i));
  // double entropy_d = 0.5 * KB * NA * diag;

  gsl_permutation * p = gsl_permutation_alloc(ndof);
  int s;
  gsl_linalg_LU_decomp(lu, p, &s);
  double lndet = gsl_linalg_LU_lndet(lu);
  double entropy = 0.5 * KB * NA * lndet;


  gsl_permutation_free(p);
  gsl_matrix_free(lu);

  return entropy;
}

double freq_etc(gsl_matrix * cov, gsl_vector * av, gsl_vector * mass, unsigned long confs, int ndof, double temperature) {
  const double dbl_confs = double(confs);
  
  // calculate factor: mu*k*T*e^2/hbar^2 (+conversion from nm to m)
  //double mukte2h2 = temperature * 1e-18 * MU * KB * E * E / (HBAR * HBAR);
  double mufac = 1e-18 * MU;
  double kT = temperature * KB;

  gsl_matrix * lu = gsl_matrix_alloc(ndof, ndof);

  double temp;
  // off-diagonal elements
  for (int i = 0; i < ndof; i++) {
    for (int j = i + 1; j < ndof; j++) {
      temp = gsl_matrix_get(cov, i, j) / dbl_confs - gsl_vector_get(av, i) * gsl_vector_get(av, j) / (dbl_confs * dbl_confs);
      temp = temp * mufac;
      // for a faster program -> use a sqrtmass_vector
      temp = temp * sqrt(gsl_vector_get(mass, i)) *
              sqrt(gsl_vector_get(mass, j));
      gsl_matrix_set(lu, i, j, temp);
      gsl_matrix_set(lu, j, i, temp);
    }
  }
  // diagonal elements
  for (int i = 0; i < ndof; i++) {
    temp = gsl_vector_get(av, i) / dbl_confs;
    temp = -temp * temp + gsl_matrix_get(cov, i, i) / dbl_confs;
    temp = temp * mufac * gsl_vector_get(mass, i);
    gsl_matrix_set(lu, i, i, temp);
  }


  gsl_vector *eval = gsl_vector_alloc(ndof);
  gsl_matrix *evec = gsl_matrix_alloc(ndof, ndof);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(ndof);
  gsl_eigen_symmv(lu, eval, evec, w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  double entropy = 0.0;
  double ho_kT;
  int nr_evals = 0;
  for (int i = 0; i < ndof; i++) {
    temp = sqrt(kT * gsl_vector_get(eval, i)) / HBAR;

    if (temp > 1.0e-10) {
      nr_evals++;
      ho_kT = 1.0 / temp;
      temp = KB * NA * (ho_kT / (exp(ho_kT) - 1.0) - log(1.0 - exp(-ho_kT)));
      entropy += temp;
    }
  }

  gsl_eigen_symmv_free(w);
  gsl_matrix_free(lu);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);

  return entropy;
}
