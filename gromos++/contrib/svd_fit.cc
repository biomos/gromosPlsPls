/**
 * @file svd_fit.cc
 * carry out an svd fit of a(n ensemble) of structures to some
 * experimental data (RDCs, 3J-values...)
 */
/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor svd_fit
 * @section svd_fit svd fit a(n ensemble) of structures to some data
 * @author @ref ja
 * @date 23.11.2009
 *
 * PROGRAM DESCRIPTION
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@type</td><td>&lt;data type: RDC, JVAL&gt; </td></tr>
 * <tr><td> \@fitspec</td><td>&lt;file containing data to fit to&gt; </td></tr>
 * <tr><td> [\@calcspec</td><td>&lt;file containing data to back-calculate (if different to @fitspec)&gt;] </td></tr>
 * <tr><td> [\@align</td><td>&lt;reference coordinates to align structures to before fitting (otherwise first frame of trj is used)&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref Atomspecifier "atoms" to consider for alignment&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;time and dt&gt;] </td></tr>
 * <tr><td> [\@timespec</td><td>&lt;timepoints to consider for the fit: ALL (default), EVERY or SPEC (if time-series)&gt;] </td></tr>
 * <tr><td> [\@timepts</td><td>&lt;timepoints to consider for the fit (if time-series and timespec EVERY or SPEC)&gt;] </td></tr>
 * <tr><td> [\@weights</td><td>&lt;file containing weights for particular frames of the trajectory&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   rdc
     @topo       ex.top
     @pbc        r
     @type       RDC
     @fitspec    ex.fit
     @atoms      1:CA
     @time       0 0.02
     @traj       ex.trj
   @endverbatim

 * <hr>
 */

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Time.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/RdcFuncs.h"
#include "../src/utils/JvalFuncs.h"
#include "../src/utils/PropertyContainer.h"


using namespace args;
using namespace bound;
using namespace fit;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace std;
using namespace utils;

// function for skipping time-points
bool use_this_frame(int i, string const & timespec, vector<int> const & timepts,
        unsigned int & timesComputed, bool & done);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "type" << "fitspec" << "calcspec" << "align" <<
         "atoms"  << "time" << "dataspec" << "timespec" << "timepts" << "weights" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>]\n";
  usage += "\t@type        <data type: RDC, JVAL>\n";
  usage += "\t@fitspec     <file containing data to fit to>\n";
  usage += "\t[@calcspec   <file containing data to backcalculate>] (if different to @fitspec)\n";
  usage += "\t[@align      <reference coordinates to align structures to before fitting>] (default: first frame of trj) (only required for RDC)\n";
  usage += "\t@atoms       <atoms to consider for alignment>\n";
  usage += "\t[@time       <t> <dt>] (optional; will only print time series if given)\n";
  usage += "\t[@timespec   <timepoints at which to compute the SASA: ALL (default), EVERY or SPEC>]\n";
  usage += "\t[@timepts    <timepoints at which to compute the SASA>] (if timespec EVERY or SPEC)\n";
  usage += "\t[@weights    <file containing weights for particular frames of the trajectory>]\n";
  usage += "\t@traj        <trajectory files>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // first check the data type
    string datatype;
    if (args.count("type") > 0) {
      datatype = args["type"];
    } else {
      throw gromos::Exception("svd_fit",
              "you must specify a data type");
    }

    // read topology and initialise systems
    InTopology it(args["topo"]);
    System sys(it.system());

    // get time
    Time time(args);

    // parse timespec
    string timespec = "ALL";
    vector<int> timepts;
    if (args.count("timespec") > 0) {
      timespec = args["timespec"];
      if (timespec != "ALL" && timespec != "EVERY" && timespec != "SPEC")
        throw gromos::Exception("svd_fit",
              "timespec format " + timespec + " unknown.\n");
      if (timespec == "EVERY" || timespec == "SPEC") {
        for (Arguments::const_iterator it = args.lower_bound("timepts");
                it != args.upper_bound("timepts"); ++it) {
          int bla = atoi(it->second.c_str());
          timepts.push_back(bla);
        }
        if (timepts.size() == 0) {
          throw gromos::Exception("svd_fit",
                  "if you give EVERY or SPEC you have to use "
                  "@timepts as well");
        }
        if (timepts.size() != 1 && timespec == "EVERY") {
          throw gromos::Exception("svd_fit",
                  "if you give EVERY you have to give exactly"
                  " one number with @timepts");
        }
      }
    }

    // only do alignment for RDCs
    if (datatype == "RDC") {
      System refSys(it.system());

      // parse boundary conditions for refSys
      Boundary *pbc = BoundaryParser::boundary(refSys, args);
      // parse gather method
      Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

      // read reference coordinates for alignment
      InG96 ic;
      if (args.count("align") > 0)
        ic.open(args["align"]);
      else
        if (args.count("traj") > 0)
        ic.open(args.lower_bound("traj")->second);
      ic >> refSys;
      ic.close();

      // check that there is a box block in the refSys
      if (refSys.hasBox == false && pbc->type() != 'v')
        throw gromos::Exception("svd_fit",
              "If pbc != v you have to give a box block "
              "in the reference system as well.");

      // and that we managed to read coordinates from it
      if (!refSys.hasPos)
        throw gromos::Exception("svd_fit",
              "Unable to read POSITION(RED) block from "
              "reference positions file.");

      // gather reference system
      (*pbc.*gathmethod)();
      delete pbc;

      Reference refalign(&refSys);

      // get alignment atoms
      AtomSpecifier alignmentatoms(refSys);
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        alignmentatoms.addSpecifier(spec);
      }
      if (alignmentatoms.size() == 0)
        throw gromos::Exception("svd_fit", "Must specify alignment atoms for RDC data!");
      refalign.addAtomSpecifier(alignmentatoms);

      // and set up some things for the alignment later
      Vec cog = PositionUtils::cog(refSys, refalign);
      RotationalFit rf(&refalign);

      // now parse boundary conditions for sys
      pbc = BoundaryParser::boundary(sys, args);

      // access RDC-specific io and fitting functions
      RdcFuncs RDCTools(sys, args);

      // check if we want to weight various frames differently
      double w = 1.0; // this is the default
      bool using_weights = false;
      vector<RDCWeights::weights> weights; // for storing weights
      int num_weights; // number of weights read from file
      if (args.count("weights") > 0) {
        using_weights = true;
        Ginstream wf(args["weights"]);
        vector<string> buffer;
        wf.getblock(buffer);

        RDCTools.read_weights(buffer, weights);
        num_weights = weights.size();

        wf.close();
      }

      // read in the fitting data from file
      Ginstream sf(args["fitspec"]);
      vector<string> buffer;
      // first search for RDCVALRESSPEC block
      bool found_rdc = false;
      while (!sf.stream().eof() && !found_rdc) {
        sf.getblock(buffer);
        if (buffer[0] == "RDCVALRESSPEC")
          found_rdc = true;
      }
      if (!found_rdc)
        throw gromos::Exception("svd_fit", "rdc file does not contain an RDCVALRESSPEC block!");

      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("svd_fit", "RDC file is corrupted. No END in "
              + buffer[0] + " block. Got\n" + buffer[buffer.size() - 1]);

      // declare some variables
      int nfit, nbc;
      vector<RDCData::rdcparam> rdc_fit;
      vector<RDCData::rdcparam> rdc_bc;

      gsl_vector *rdc_vec, *S, *Stmp, *work;
      gsl_matrix *coef_mat, *A, *V;

      // read in the RDC data
      RDCTools.read_rdc(buffer, sys, rdc_fit);
      nfit = rdc_fit.size();
      sf.close();
      // compute rij and so dmax from reference system (constant if bonded)
      RDCTools.calc_dmax(refSys, rdc_fit, nfit);

      // allocate sizes
      coef_mat = gsl_matrix_alloc(nfit, 5);
      gsl_matrix_set_zero(coef_mat);

      A = gsl_matrix_alloc(nfit, 5);
      V = gsl_matrix_alloc(5, 5);
      S = gsl_vector_alloc(5);
      Stmp = gsl_vector_alloc(5);
      work = gsl_vector_alloc(5);
      rdc_vec = gsl_vector_alloc(nfit);

      // now we can fill rdc_vec
      RDCTools.fill_rdcvec(rdc_fit, rdc_vec, nfit);

      // declare arrays for back-calculated RDC data
      gsl_vector *rdc_vec_bc, *bc, *rdc_vec_bc_unnorm, *bc_unnorm;
      gsl_matrix *coef_mat_bc;

      // check if we have different data to back-calculate
      bool backcalc = false;
      if (args.count("calcspec") > 0) {
        backcalc = true;
        Ginstream sf(args["calcspec"]);
        sf.getblock(buffer);

        RDCTools.read_rdc(buffer, sys, rdc_bc);
        nbc = rdc_bc.size();

        sf.close();
      }
      // otherwise copy fit data into bc data
      else {
        rdc_bc = rdc_fit;
        nbc = rdc_bc.size();
      }

      // initialise arrays for back-calculated RDC data
      rdc_vec_bc_unnorm = gsl_vector_alloc(nfit);
      bc = gsl_vector_alloc(nbc);
      bc_unnorm = gsl_vector_alloc(nfit);
      gsl_vector_set_zero(bc);

      // start at -1 to get times right
      int num_frames = -1;
      // number of time-points for which coefficients for svd fit have been calculated
      unsigned int times_computed = 0;
      // for SPEC: so that we stop trying when all requested timepoints are written
      bool done = false;
      // number of frame weights used so far (in case number of weights != number of frames)
      int weights_used = 0;
      // sum of weights (used to normalise the coefficients)
      double weight_sum = 0;

      // loop over all trajectories
      for (Arguments::const_iterator
        iter = args.lower_bound("traj"),
              to = args.upper_bound("traj");
              iter != to; ++iter) {

        // open file
        ic.open((iter->second).c_str());
        ic.select("SOLUTE");

        // loop over single trajectory
        while (!ic.eof()) {

          ic >> sys >> time;
          if (!sys.hasPos)
            throw gromos::Exception("rmsd",
                  "Unable to read POSITION(RED) block from "
                  "trajectory file.");

          // gather
          (*pbc.*gathmethod)();

          // check whether to skip this frame or not
          num_frames++;
          if (use_this_frame(num_frames, timespec, timepts, times_computed, done)) {

            // align to reference structure
            rf.fit(&sys);

            // get weight from file (otherwise already set to 1.0)
            if (using_weights && weights_used <= num_weights &&
                    weights[weights_used].frame == num_frames) {
              w = weights[weights_used].weight;
              weights_used++;
            }
            weight_sum += w;

            // calculate the coefficients
            RDCTools.calc_coef(sys, rdc_fit, coef_mat, nfit, w);

          }//end if use_this_frame
          if (done)
            break;

        }//end frame
      }// end io

      // scale coefficients by number of frames (sum of weights as each contribution
      // was multiplied by its weight)
      gsl_matrix_scale(coef_mat, 1. / weight_sum);

      // check if we want to back-calculate the same data that we fit to
      if (!backcalc) {
        // back-calculating same data as fitted to so just copy coefficients
        coef_mat_bc = coef_mat;
        rdc_vec_bc = rdc_vec;
      } else {
        // get bc rdc data and initialise coefficient matrix
        rdc_vec_bc = gsl_vector_alloc(nfit);
        RDCTools.fill_rdcvec(rdc_bc, rdc_vec_bc, nbc);
        coef_mat_bc = gsl_matrix_alloc(nfit, 5);
        gsl_matrix_set_zero(coef_mat_bc);
      }

      // make a copy of coef_mat because it is modified by the SVD functions
      gsl_matrix_memcpy(A, coef_mat);

      // calculate and print alignment tensor (fitting to fit data)      
      gsl_linalg_SV_decomp(A, V, Stmp, work);
      gsl_linalg_SV_solve(A, V, Stmp, rdc_vec, S);

      double Sxx = gsl_vector_get(S, 0);
      double Syy = gsl_vector_get(S, 1);
      double Szz = -Sxx - Syy;
      double Sxy = gsl_vector_get(S, 2);
      double Sxz = gsl_vector_get(S, 3);
      double Syz = gsl_vector_get(S, 4);

      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << "# Alignment tensor" << endl;
      fprintf(stdout,"# Sxx = %12.5e\n",Sxx);
      fprintf(stdout,"# Syy = %12.5e\n",Syy);
      fprintf(stdout,"# Szz = %12.5e\n",Szz);
      fprintf(stdout,"# Sxy = %12.5e\n",Sxy);
      fprintf(stdout,"# Sxz = %12.5e\n",Sxz);
      fprintf(stdout,"# Syz = %12.5e\n",Syz);

      // calculate and print Q_norm (coef_mat_bc is still empty...bc gets filled?)
      gsl_blas_dgemv(CblasNoTrans, 1.0, coef_mat_bc, S, 0., bc);
      double Q_norm = RDCTools.calc_Q(bc, rdc_vec_bc);
      fprintf(stdout,"# Normalized Q = %12.6f\n",Q_norm);

      // calculate and print Q_raw
      RDCTools.unnorm_rdcvec(rdc_bc, bc, bc_unnorm);
      RDCTools.unnorm_rdcvec(rdc_bc, rdc_vec_bc, rdc_vec_bc_unnorm);
      double Q_raw = RDCTools.calc_Q(bc_unnorm, rdc_vec_bc_unnorm);
      fprintf(stdout, "# Raw Q =        %12.6f\n", Q_raw);

      cout.setf(ios::fixed, ios::floatfield);
      cout.precision(5);

      // title for actual RDCs
      cout << setw(9) << "#   atom1" << setw(9) << "atom2"
      << setw(11) << "D_calc" << setw(11) << "D_exp"
      << setw(11) << "|D_calc|" << setw(11) << "|D_exp|"  << endl;

      // back-conversion factor to scale from ps-1 into Hz (s-1)
      double rdc_bcf = 1.0 / pico;

      // now print actual (back-calculated) RDCs
      for (int i = 0; i < nbc; i++) {

              cout << setw(5) << rdc_bc[i].i + 1
              << setw(4) << sys.mol(rdc_bc[i].mol).topology().atom(rdc_bc[i].i).name()
              << setw(5) << rdc_bc[i].j + 1
              << setw(4) << sys.mol(rdc_bc[i].mol).topology().atom(rdc_bc[i].j).name()
              << setw(11) << gsl_vector_get(bc_unnorm, i) * rdc_bcf
              << setw(11) << gsl_vector_get(rdc_vec_bc_unnorm, i) * rdc_bcf
              << setw(11) << gsl_vector_get(bc, i) * rdc_bc[i].dmax * rdc_bcf
              << setw(11) << gsl_vector_get(rdc_vec_bc, i) * rdc_bc[i].dmax * rdc_bcf
              << endl;
      }

      // clean up??
      gsl_vector_free(rdc_vec);
      gsl_matrix_free(coef_mat);


    }

    // start JVAL-specific loop
    else if (datatype == "JVAL") {

      // parse boundary conditions for sys
      Boundary *pbc = BoundaryParser::boundary(sys, args);
      // parse gather method
      Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

      // access JVAL-specific io and fitting functions
      JvalFuncs JvalTools(sys, args);

      // check if we want to weight various frames differently
      double w = 1.0; // this is the default
      bool using_weights = false;
      vector<JvalWeights::weights> weights; // for storing weights
      int num_weights; // number of weights read from file
      if (args.count("weights") > 0) {
        using_weights = true;
        Ginstream wf(args["weights"]);
        vector<string> buffer;
        wf.getblock(buffer);

        JvalTools.read_weights(buffer, weights);
        num_weights = weights.size();

        wf.close();
      }

      // read in the fitting data from file
      Ginstream sf(args["fitspec"]);
      vector<string> buffer;
      // first search for JVALRESSPEC block
      bool found_jval = false;
      while (!sf.stream().eof() && !found_jval) {
        sf.getblock(buffer);
        if (buffer[0] == "JVALRESSPEC")
          found_jval = true;
      }
      if (!found_jval)
        throw gromos::Exception("svd_fit", "JVAL file does not contain a JVALRESSPEC block!");

      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("svd_fit", "JVAL file is corrupted. No END in "
              + buffer[0] + " block. Got\n" + buffer[buffer.size() - 1]);

      // declare some variables
      int nfit, nbc;
      vector<JvalData::jvalues> jval_fit;
      vector<JvalData::jvalues> jval_bc;

      gsl_vector *jval_vec, *K, *Ktmp, *work;
      gsl_matrix *coef_mat, *A, *V;

      // read in the Jvalue data
      PropertyContainer fit_props(sys, pbc);
      double delta;
      JvalTools.read_jval(buffer, sys, jval_fit, fit_props, delta);
      nfit = jval_fit.size();
      sf.close();

      // allocate sizes
      coef_mat = gsl_matrix_alloc(nfit, 3);
      gsl_matrix_set_zero(coef_mat);

      A = gsl_matrix_alloc(nfit, 3);
      V = gsl_matrix_alloc(3, 3);
      K = gsl_vector_alloc(3);
      Ktmp = gsl_vector_alloc(3);
      work = gsl_vector_alloc(3);
      jval_vec = gsl_vector_alloc(nfit);

      // now we can fill jval_vec
      JvalTools.fill_jvalvec(jval_fit, jval_vec, nfit);

      // declare arrays for back-calculated RDC data
      gsl_vector *jval_vec_bc, *bc;
      gsl_matrix *coef_mat_bc;

      // check if we have different data to back-calculate
      bool backcalc = false;
      if (args.count("calcspec") > 0) {
        backcalc = true;
        Ginstream sf(args["calcspec"]);
        sf.getblock(buffer);

        PropertyContainer bc_props(sys, pbc);
        JvalTools.read_jval(buffer, sys, jval_bc, bc_props, delta);
        nbc = jval_bc.size();

        sf.close();
      }
      // otherwise copy fit data into bc data
      else {
        jval_bc = jval_fit;
        nbc = jval_bc.size();
      }

      // initialise arrays for back-calculated JVAL data
      bc = gsl_vector_alloc(nbc);
      gsl_vector_set_zero(bc);

      // start at -1 to get times right
      int num_frames = -1;
      // number of time-points for which coefficients for svd fit have been calculated
      unsigned int times_computed = 0;
      // for SPEC: so that we stop trying when all requested timepoints are written
      bool done = false;
      // number of frame weights used so far (in case number of weights != number of frames)
      int weights_used = 0;
      // sum of weights (used for normalisation)
      double weight_sum = 0;

      // loop over all trajectories
      for (Arguments::const_iterator
        iter = args.lower_bound("traj"),
              to = args.upper_bound("traj");
              iter != to; ++iter) {

        // open file
        InG96 ic;
        ic.open((iter->second).c_str());
        ic.select("SOLUTE");

        // loop over single trajectory
        while (!ic.eof()) {

          ic >> sys >> time;
          if (!sys.hasPos)
            throw gromos::Exception("svd_fit",
                  "Unable to read POSITION(RED) block from "
                  "trajectory file.");

          // gather
          (*pbc.*gathmethod)();

          // check whether to skip this frame or not
          num_frames++;
          if (use_this_frame(num_frames, timespec, timepts, times_computed, done)) {

            // get weight from file (otherwise already set to 1.0)
            if (using_weights && weights_used <= num_weights &&
                    weights[weights_used].frame == num_frames) {
              w = weights[weights_used].weight;
              weights_used++;
            }
            // w is alredy 0.0 (if not reset above)
            weight_sum += w;

            // calculate the coefficients
            JvalTools.calc_coef(sys, fit_props, coef_mat, nfit, w, delta);

          }//end if use_this_frame
          if (done)
            break;

        }//end frame
      }// end io

      // scale coefficients by number of frames (sum of weights as each contribution
      // was multiplied by its weight)
      gsl_matrix_scale(coef_mat, 1. / weight_sum);

      // check if we want to back-calculate the same data that we fit to
      if (!backcalc) {
        // back-calculating same data as fitted to so just copy coefficients
        coef_mat_bc = coef_mat;
        jval_vec_bc = jval_vec;
      } else {
        // get bc jval data and initialise coefficient matrix
        jval_vec_bc = gsl_vector_alloc(nfit);
        JvalTools.fill_jvalvec(jval_bc, jval_vec_bc, nbc);
        coef_mat_bc = gsl_matrix_alloc(nfit, 3);
        gsl_matrix_set_zero(coef_mat_bc);
      }

      // make a copy of coef_mat because it is modified by the SVD functions
      gsl_matrix_memcpy(A, coef_mat);

      // calculate and print Karplus parameters (fitting to fit data)
      gsl_linalg_SV_decomp(A, V, Ktmp, work);
      gsl_linalg_SV_solve(A, V, Ktmp, jval_vec, K);

      // back-conversion factor to scale from ps-1 into Hz (s-1): pointless
      //double bcf = 1.0 / pico;
      //double K_A = gsl_vector_get(K, 0) * bcf;
      //double K_B = gsl_vector_get(K, 1) * bcf;
      //double K_C = gsl_vector_get(K, 2) * bcf;

      double K_A = gsl_vector_get(K, 0);
      double K_B = gsl_vector_get(K, 1);
      double K_C = gsl_vector_get(K, 2);

      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << "# Optimised Karplus parameters:" << endl;
      fprintf(stdout,"# A = %12.5e\n",K_A);
      fprintf(stdout,"# B = %12.5e\n",K_B);
      fprintf(stdout,"# C = %12.5e\n",K_C);

      // back-calculate J-values (this fills bc with coef_mat_bc * K)
      gsl_blas_dgemv(CblasNoTrans, 1.0, coef_mat_bc, K, 0., bc);
      // calculate and print Q
      double Q = JvalTools.calc_Q(bc, jval_vec_bc);
      cout << "# Goodness of fit:" << endl;
      fprintf(stdout,"# Q = %12.6f\n",Q);

      cout.setf(ios::fixed, ios::floatfield);
      cout.precision(5);

      // title for back-calculated Jvalues
      cout << setw(5) << "#   i" << setw(5) << "j" << setw(5) << "k" << setw(5) << "l"
      << setw(11) << "J_calc" << setw(11) << "J_exp" << endl;

      // now print actual (back-calculated) Jvalues
      for (int i = 0; i < nbc; i++) {

        cout << setw(5) << jval_bc[i].i + 1
                << setw(5) << jval_bc[i].j + 1
                << setw(5) << jval_bc[i].k + 1
                << setw(5) << jval_bc[i].l + 1
                //<< setw(11) << gsl_vector_get(bc, i) * bcf
                //<< setw(11) << gsl_vector_get(jval_vec_bc, i) * bcf
                << setw(11) << gsl_vector_get(bc, i)
                << setw(11) << gsl_vector_get(jval_vec_bc, i)
                << endl;
      }

      // clean up??
      gsl_vector_free(jval_vec);
      gsl_matrix_free(coef_mat);

    }

    else {
      throw gromos::Exception("svd_fit", "Data types other than RDC "
              "not implemented yet");
    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

// function to decide whether to use this frame or not
bool use_this_frame(int i, std::string const & timespec, vector<int> const & timepts,
        unsigned int & timesComputed, bool & done) {
  if (timespec == "ALL") {
    ++timesComputed;
    return true;
  } else if (timespec == "EVERY" && i % timepts[0] == 0) {
    ++timesComputed;
    return true;
  } else if (timespec == "SPEC") {
    for (unsigned int j = 0; j < timepts.size(); ++j) {
      if (timepts[j] == i) {
        ++timesComputed;
        if (timesComputed == timepts.size())
          done = true;
        return true;
      } // compute
    } // times
  }
  return false;
}

