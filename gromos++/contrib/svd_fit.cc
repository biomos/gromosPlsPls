/**
 * @file svd_fit.cc
 * carry out an svd fit of a(n ensemble of) structures to some
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
 * <tr><td> \@fitrij</td><td>&lt;inter-nuclear distances to use for fit data: SPEC read from @fitspec, INIT calculate from @align, ALL calculate from each structure in @traj&gt; </td></tr>
 * <tr><td> [\@bcspec</td><td>&lt;file containing data to back-calculate (if different to @fitspec)&gt;] </td></tr>
 * <tr><td> [\@bcrij</td><td>&lt;inter-nuclear distances to use for bc data: SPEC, INIT or ALL]&gt; </td></tr>
 * <tr><td> [\@scale</td><td>&lt;scale non-NH RDCs to match NH RDCs using gyroN, gyroH (*10^6 rad/Ts = Hz/s) and rNH (nm) (default: -27.12, 267.52, 0.1)] (so that all RDC types are weighted similarly when fitting to multiple types simultaneously)&gt; </td></tr>
 * <tr><td> [\@align</td><td>&lt;reference coordinates to align structures to before fitting (otherwise first frame of trj is used)&gt;] </td></tr>
 * <tr><td> [\@atoms</td><td>&lt;@ref Atomspecifier "atoms" to consider for alignment (compulsory if @type is RDC)]&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;time and dt&gt;] </td></tr>
 * <tr><td> [\@timespec</td><td>&lt;timepoints to consider for the fit: ALL (default), EVERY or SPEC (if time-series)&gt;] </td></tr>
 * <tr><td> [\@timepts</td><td>&lt;timepoints to consider for the fit (if time-series and timespec EVERY or SPEC)&gt;] </td></tr>
 * <tr><td> [\@weights</td><td>&lt;file containing weights for particular frames of the trajectory&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
 svd_fit
     @topo       ex.top
     @pbc        r
     @type       RDC
     @fitspec    ex.fit
     @fitrij     ALL
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
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
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

// -- function for skipping time-points
bool use_this_frame(int i, string const & timespec, vector<int> const & timepts,
        unsigned int & timesComputed, bool & done);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "type" << "fitspec" << "fitrij" 
         << "bcspec" << "bcrij" << "scale" << "align" << "atoms"  << "time"
         << "dataspec" << "timespec" << "timepts" << "weights" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>]\n";
  usage += "\t@type        <data type: RDC, JVAL>\n";
  usage += "\t@fitspec     <file containing data to fit to>]\n";
  usage += "\t[@fitrij     <inter-nuclear distances to use for fit data: SPEC read from @fitspec, INIT calculate from @align, ALL calculate from each structure in @traj>] (only for @type RDC)\n";
  usage += "\t[@bcspec     <file containing data to backcalculate>] (if different to @fitspec)\n";
  usage += "\t[@bcrij      <inter-nuclear distances to use for bc data: SPEC, INIT or ALL>]\n";
  usage += "\t[@scale      <scale non-NH RDCs to match NH RDCs using gyroN, gyroH (*10^6 rad/Ts = Hz/s) and rNH (nm) (default: 27.12, 267.52, 0.1)>] (so that all RDC types are weighted similarly when fitting to multiple types simultaneously)\n";
  usage += "\t[@align      <reference coordinates to align structures to before fitting>] (default: first frame of trj) (only required for RDC)\n";
  usage += "\t[@atoms      <atoms to consider for alignment>] (only for @type RDC)\n";
  usage += "\t[@time       <t> <dt>] (optional; will only print time series if given)\n";
  usage += "\t[@timespec   <timepoints to consider for the fit: ALL (default), EVERY or SPEC>]\n";
  usage += "\t[@timepts    <timepoints to consider for the fit>] (if timespec EVERY or SPEC)\n";
  usage += "\t[@weights    <file containing weights for particular frames of the trajectory>]\n";
  usage += "\t@traj        <trajectory files>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // -- read topology and initialise system
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    // -- get time
    Time time(args);

    // -- parse timespec
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

    // -- check data type
    string datatype;
    if (args.count("type") > 0) {
      datatype = args["type"];
      //args.getValue<std::string>("type", true);
    } else {
      throw gromos::Exception("svd_fit",
              "you must specify a data type");
    }

    // -- svd_fit for RDCs
    if (datatype == "RDC") {

      // declare existence of reference system
      System refSys(it.system());
      // parse boundary conditions for refSys
      Boundary *pbc = BoundaryParser::boundary(refSys, args);
      // parse gather method
      Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

      // -- read reference coordinates for alignment
      InG96 ic;
      if (args.count("align") > 0) {
        ic.open(args["align"]);
      } else if (args.count("traj") > 0) {
        ic.open(args.lower_bound("traj")->second);
      } else {
        throw gromos::Exception("svd_fit",
                "You need to give some coordinates!");
      }
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
      // and declare it as a reference, with new name refalign
      Reference refalign(&refSys);

      // -- get alignment atoms
      AtomSpecifier alignmentatoms(refSys);
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        alignmentatoms.addSpecifier(spec);
      }
      if (alignmentatoms.size() == 0) {
        throw gromos::Exception("svd_fit", "Must specify alignment atoms for RDC data!");
      } else if (alignmentatoms.size() < 4) {
        throw gromos::Exception("svd_fit", "At least 4 alignment atoms required for rotational fit\n");
      }
      refalign.addAtomSpecifier(alignmentatoms);
      RotationalFit rf(&refalign); // rmsd

      // shift coordinates of reference so that CoM is at origin
      PositionUtils::shiftToCom(&refSys);
      // now parse boundary conditions for sys
      pbc = BoundaryParser::boundary(sys, args);
      // access RDC-specific io and fitting functions
      RdcFuncs RDCTools(sys, args);

      // -- check if we want to weight various frames differently
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

      // -- read in the fitting data from file
      Ginstream sf(args["fitspec"]);
      vector<string> buffer;
      // first search for RDCRESSPEC block
      bool found_rdc = false;
      while (!sf.stream().eof() && !found_rdc) {
        sf.getblock(buffer);
        if (buffer[0] == "RDCRESSPEC")
          found_rdc = true;
      }
      if (!found_rdc)
        throw gromos::Exception("svd_fit", "RDC fit file does not contain an RDCRESSPEC block!");
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("svd_fit", "RDC fit file is corrupted. No END in "
              + buffer[0] + " block. Got\n" + buffer[buffer.size() - 1]);

      // check how to get the inter-nuclear distances
      string fit_rij;
      if (args.count("fitrij") > 0) {
        fit_rij = args["fitrij"];
        if (fit_rij != "SPEC" && fit_rij != "INIT" && fit_rij != "ALL")
          throw gromos::Exception("svd_fit",
                "fitrij format " + fit_rij + " unknown.\n");
      }

      // do we scale nonNH RDCs?
      bool scale = false;
      if (args.count("scale") >= 0) {
        scale = true;
      }
      vector<double> scalearg = args.getValues<double>("scale", 3, false,
          Arguments::Default<double>() << 27.12 << 267.52 << 0.1);
      double gyroN = scalearg[0];
      double gyroH = scalearg[1];
      double rNH = scalearg[2];

      // convert into gromos units: hard-coded because gromos++ uses slightly different values to MD++
      // GROMOS++
      //const double atomic_mass_unit = gmath::physConst.get_atomic_mass_unit();
      //const double elementary_charge = gmath::physConst.get_elementary_charge();
      // MD++ (from Wilfred's notes ~1970)
      const double atomic_mass_unit = 1.6605655e-27;
      const double elementary_charge = 1.6021892e-19;
      const double g_cf = (atomic_mass_unit / elementary_charge) * 1.0e6;
      gyroN *= g_cf;
      gyroH *= g_cf;

      // container for fit data
      unsigned int nfit;
      vector<RDCData::rdcparam> fit_dat;

      // read in the RDC data for fitting
      RDCTools.read_rdc(buffer, sys, fit_dat, fit_rij, true);
      nfit = fit_dat.size();
      sf.close();

      // compute rij if INIT
      if (fit_rij == "INIT") {
        RDCTools.init_rij(refSys, fit_dat);
      }

      // check if we have different data to back-calculate
      bool backcalc = false;
      string bc_rij;
      unsigned int nbc = nfit;
      vector<RDCData::rdcparam> bc_dat;
      gsl_vector *exp_bc_norm;
      if (args.count("bcspec") > 0) {
        backcalc = true;
        // check how to get the inter-nuclear distances
        if (args.count("bcrij") > 0) {
          bc_rij = args["bcrij"];
          if (bc_rij != "SPEC" && bc_rij != "INIT" && bc_rij != "ALL")
            throw gromos::Exception("svd_fit",
                  "bcrij format " + bc_rij + " unknown.\n");
        }

        // read in the data to back-calculate from file
        Ginstream sf(args["bcspec"]);
        vector<string> buffer;
        // first search for RDCRESSPEC block
        bool found_rdc = false;
        while (!sf.stream().eof() && !found_rdc) {
          sf.getblock(buffer);
          if (buffer[0] == "RDCRESSPEC")
            found_rdc = true;
        }
        if (!found_rdc)
          throw gromos::Exception("svd_fit", "RDC bc file does not contain an RDCRESSPEC block!");
        if (buffer[buffer.size() - 1].find("END") != 0)
          throw gromos::Exception("svd_fit", "RDC bc file is corrupted. No END in "
                + buffer[0] + " block. Got\n" + buffer[buffer.size() - 1]);

        RDCTools.read_rdc(buffer, sys, bc_dat, bc_rij, false);
        nbc = bc_dat.size();
        sf.close();

        // compute rij from initial structure (if INIT)
        if (bc_rij == "INIT") {
          RDCTools.init_rij(refSys, bc_dat);
        }
      }

      // start at -1 to get times right. this is for the timespec.
      int num_frames = -1;
      // and this one is for normalising rij
      double num_rij = 0.0;
      // number of time-points for which coefficients for svd fit have been calculated
      unsigned int times_computed = 0;
      // for SPEC: so that we stop trying when all requested timepoints are written
      bool done = false;
      // number of frame weights used so far (in case number of weights != number of frames)
      int weights_used = 0;
      // sum of weights (used to normalise the coefficients)
      double weight_sum = 0;
      // gsl used within trajectory loop
      gsl_matrix *coef_mat = gsl_matrix_calloc(nfit, 5);
      // two sets of coefficients for back-calculation in case of side-chain NH RDCs
      gsl_matrix *coef_mat_bc_j = gsl_matrix_calloc(nbc, 5);
      gsl_matrix *coef_mat_bc_k = gsl_matrix_calloc(nbc, 5);

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
            throw gromos::Exception("svd_fit",
                  "Unable to read POSITION(RED) block from trajectory file.");
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

            // calculate the coefficients for the RDCs to fit to
            // and the inter-nuclear distances if fit_rij == ALL
            // and find HA coords if there are CA-HA RDCs
            RDCTools.calc_coef_fit(sys, fit_dat, coef_mat, nfit, w, fit_rij);

            // and for the RDCs to back-calculate (if different)
            if (backcalc) {
              RDCTools.calc_coef_bc(sys, bc_dat, coef_mat_bc_j, coef_mat_bc_k,
                      nbc, w, bc_rij);
            }

            // update num_rij
            num_rij += 1.0;

          }//end if use_this_frame
          if (done)
            break;

        }//end frame

      }// end io

      // scale coefficients by number of frames (sum of weights as each contribution
      // was multiplied by its weight)
      gsl_matrix_scale(coef_mat, 1. / weight_sum);

      // scale rij by number of frames (if fit_rij == ALL)
      if (fit_rij == "ALL") {
        RDCTools.scale_rij(fit_dat, num_rij);
      }
      // compute Dmax (includes scaling of prefactor for nonNH RDCs if requested)
      RDCTools.calc_dmax8(fit_dat, sys, scale, gyroN, gyroH, rNH);
      // normalise the RDCs, and put them in a gsl vector
      gsl_vector *exp_fit_norm;
      exp_fit_norm = gsl_vector_alloc(nfit);
      RDCTools.fill_rdcvec_norm(fit_dat, exp_fit_norm);

      // if we are back-calculating the same data that we fit to, we can just
      // copy the coefficients and RDCs across
      exp_bc_norm = gsl_vector_alloc(nbc);
      if (!backcalc) {
        bc_dat = fit_dat;
        coef_mat_bc_j = coef_mat;
        // (coef_mat_bc_k can stay all zeros as we can't have fitted to side-chain NH)
        exp_bc_norm = exp_fit_norm;
      }
      // if not, the coefficients were already calculated
      else {
        // scale the coefficients
        gsl_matrix_scale(coef_mat_bc_j, 1. / weight_sum);
        gsl_matrix_scale(coef_mat_bc_k, 1. / weight_sum);
        // scale rij (if ALL)
        if (bc_rij == "ALL") {
          RDCTools.scale_rij(bc_dat, num_rij);
        }
        // compute Dmax (includes scaling of prefactor for nonNH RDCs if requested)
        RDCTools.calc_dmax8(bc_dat, sys, scale, gyroN, gyroH, rNH);
        // fill exp_bc_sc with normalised RDCs
        RDCTools.fill_rdcvec_norm(bc_dat, exp_bc_norm);
      }

      // -- alignment tensor --

      // set up gsl for svd fitting
      gsl_matrix *A = gsl_matrix_alloc(nfit, 5);
      gsl_matrix *V = gsl_matrix_alloc(5, 5);
      gsl_vector *S = gsl_vector_alloc(5);
      gsl_vector *Stmp = gsl_vector_alloc(5);
      gsl_vector *swork = gsl_vector_alloc(5);
      // make a copy of coef_mat because it is modified by the SVD functions
      gsl_matrix_memcpy(A, coef_mat);
      // calculate alignment tensor (fitting to fit data)      
      gsl_linalg_SV_decomp(A, V, Stmp, swork);
      gsl_vector_free(swork);
      gsl_linalg_SV_solve(A, V, Stmp, exp_fit_norm, S);
      gsl_vector_free(Stmp);
      gsl_matrix_free(A);
      gsl_matrix_free(V);
      // extract components
      double Sxx = gsl_vector_get(S, 0);
      double Syy = gsl_vector_get(S, 1);
      double Szz = -Sxx - Syy;
      double Sxy = gsl_vector_get(S, 2);
      double Sxz = gsl_vector_get(S, 3);
      double Syz = gsl_vector_get(S, 4);
      // print alignment tensor
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << "# Alignment tensor" << endl;
      fprintf(stdout,"# Sxx = %12.5e\n",Sxx);
      fprintf(stdout,"# Syy = %12.5e\n",Syy);
      fprintf(stdout,"# Szz = %12.5e\n",Szz);
      fprintf(stdout,"# Sxy = %12.5e\n",Sxy);
      fprintf(stdout,"# Sxz = %12.5e\n",Sxz);
      fprintf(stdout,"# Syz = %12.5e\n",Syz);

      // -- eigenvectors and eigenvalues of S --
      // (i.e. diagonalise the alignment tensor, S)

      // put S into matrix format (rows, cols):
      gsl_matrix *Smat = gsl_matrix_alloc(3,3);
      gsl_matrix_set(Smat,0,0,Sxx);
      gsl_matrix_set(Smat,0,1,Sxy);
      gsl_matrix_set(Smat,0,2,Sxz);
      gsl_matrix_set(Smat,1,0,Sxy); // == Syx
      gsl_matrix_set(Smat,1,1,Syy);
      gsl_matrix_set(Smat,1,2,Syz);
      gsl_matrix_set(Smat,2,0,Sxz); // == Szx
      gsl_matrix_set(Smat,2,1,Syz); // == Szy
      gsl_matrix_set(Smat,2,2,Szz);      
      // allocate space (O4n, where n is the size given (for nxn matrix))
      gsl_eigen_symmv_workspace *ework =  gsl_eigen_symmv_alloc(3);
      // set up eigenvalue vector and eigenvector matrix
      gsl_matrix *evec = gsl_matrix_alloc(3,3);
      gsl_vector *eval = gsl_vector_alloc(3);
      // compute the eigenvalues and eigenvectors
      gsl_eigen_symmv(Smat, eval, evec, ework);
      // free space
      gsl_eigen_symmv_free(ework);
      gsl_matrix_free(Smat);
      // sort the eigenvalues in descending order of magnitude
      gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);
      // get the eigenvalues
      double eval_zz = gsl_vector_get(eval, 0);
      double eval_yy = gsl_vector_get(eval, 1);
      double eval_xx = gsl_vector_get(eval, 2);
      // print the first three eigenvalues
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << "# Eigenvalues" << endl;
      fprintf(stdout, "# Sxx_d = %12.5e\n", eval_xx);
      fprintf(stdout, "# Syy_d = %12.5e\n", eval_yy);
      fprintf(stdout, "# Szz_d = %12.5e\n", eval_zz);
      // get the eigenvectors (stored in the columns of evec) (rows, cols)
      double evec_zx = gsl_matrix_get(evec, 0, 0);
      double evec_zy = gsl_matrix_get(evec, 1, 0);
      double evec_zz = gsl_matrix_get(evec, 2, 0);
      double evec_yx = gsl_matrix_get(evec, 0, 1);
      double evec_yy = gsl_matrix_get(evec, 1, 1);
      double evec_yz = gsl_matrix_get(evec, 2, 1);
      double evec_xx = gsl_matrix_get(evec, 0, 2);
      double evec_xy = gsl_matrix_get(evec, 1, 2);
      double evec_xz = gsl_matrix_get(evec, 2, 2);
      // free space
      gsl_vector_free(eval);
      gsl_matrix_free(evec);
      // print the first three eigenvectors
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << "# Eigenvectors   x_coor   y_coor   z_coor" << endl;
      fprintf(stdout, "# x-axis %12.5e %12.5e %12.5e\n", evec_xx, evec_xy, evec_xz);
      fprintf(stdout, "# y-axis %12.5e %12.5e %12.5e\n", evec_yx, evec_yy, evec_yz);
      fprintf(stdout, "# z-axis %12.5e %12.5e %12.5e\n", evec_zx, evec_zy, evec_zz);

      // -- compute alignment parameters --
      double Da = 0.5 * eval_zz;
      double xydiff = eval_xx - eval_yy;
      double Dr = (xydiff) / 3.0;
      double Aa = eval_zz;
      double Ar = (2.0 / 3.0) * xydiff;
      double R = Ar / Aa; // rhombicity
      // print alignment parameters
      cout << "# Alignment parameters" << endl;
      fprintf(stdout, "# Da = %12.5e\n", Da);
      fprintf(stdout, "# Dr = %12.5e\n", Dr);
      fprintf(stdout, "# Aa = %12.5e\n", Aa);
      fprintf(stdout, "# Ar = %12.5e\n", Ar);
      fprintf(stdout, "# R = %12.5e (rhombicity)\n", R);

      // -- compute Euler angles --
      Vec xeigen(evec_xx, evec_xy, evec_xz);
      Vec yeigen(evec_yx, evec_yy, evec_yz);
      Vec zeigen(evec_zx, evec_zy, evec_zz);
      // each vector is a column of the matrix
      Matrix eigen(xeigen, yeigen, zeigen); // all we need, as "axes" matrix is identity matrix
      // x-convention Euler angles
      double alpha1, alpha2, beta1, beta2, gamma1, gamma2;
      RDCTools.euler_x(eigen, alpha1, alpha2, beta1, beta2, gamma1, gamma2);
      cout << "# Euler angles for rotation of molecule frame into PAS of alignment tensor\n"
              "# according to \"x-convention\" (anti-clockwise rotation around z,x',z'')" << endl;
      cout << "# alpha (around z-axis)   1: " << alpha1 << "\t2: " << alpha2 << endl;
      cout << "# beta  (around x'-axis)  1: " << beta1 << "\t2: " << beta2 << endl;
      cout << "# gamma (around z''-axis) 1: " << gamma1 << "\t2: " << gamma2 << endl;
      // y-convention Euler angles
      RDCTools.euler_y(eigen, alpha1, alpha2, beta1, beta2, gamma1, gamma2);
      cout << "# Euler angles for rotation of molecule frame into PAS of alignment tensor\n"
              "# according to \"y-convention\" (anti-clockwise rotation around z,y',z'')" << endl;
      cout << "# alpha (around z-axis)   1: " << alpha1 << "\t2: " << alpha2 << endl;
      cout << "# beta  (around y'-axis)  1: " << beta1 << "\t2: " << beta2 << endl;
      cout << "# gamma (around z''-axis) 1: " << gamma1 << "\t2: " << gamma2 << endl;
      // pitch-roll-yaw Euler angles
      double psi1, psi2, theta1, theta2, phi1, phi2;
      RDCTools.euler_pry(eigen, psi1, psi2, theta1, theta2, phi1, phi2);
      cout << "# Euler angles for rotation of molecule frame into PAS of alignment tensor\n"
              "# according to pitch-roll-yaw convention (anti-clockwise rotation around x,y,z)" << endl;
      cout << "# psi   (around x-axis)   1: " << psi1 << "\t2: " << psi2 << endl;
      cout << "# theta (around y-axis)   1: " << theta1 << "\t2: " << theta2 << endl;
      cout << "# phi   (around z-axis)   1: " << phi1 << "\t2: " << phi2 << endl;

      // -- back-calculate the bc-RDCs using S
      gsl_vector *calc_bc_norm_j = gsl_vector_calloc(nbc);
      gsl_vector *calc_bc_norm_k = gsl_vector_calloc(nbc);
      gsl_blas_dgemv(CblasNoTrans, 1.0, coef_mat_bc_j, S, 0., calc_bc_norm_j);
      gsl_blas_dgemv(CblasNoTrans, 1.0, coef_mat_bc_k, S, 0., calc_bc_norm_k);

      // -- sum the back-calculated RDCs
      gsl_vector *calc_bc_norm = gsl_vector_alloc(nbc);
      RDCTools.sum_rdc(nbc, calc_bc_norm_j, calc_bc_norm_k, calc_bc_norm);
      
      // -- un-normalise the RDCs (in gsl format)
      gsl_vector *calc_bc = gsl_vector_alloc(nbc);
      gsl_vector *exp_bc = gsl_vector_alloc(nbc);
      RDCTools.unnorm_rdc(bc_dat, calc_bc_norm, calc_bc);
      RDCTools.unnorm_rdc(bc_dat, exp_bc_norm, exp_bc);

      // -- calculate and print Q-values, R-values and rmsd
      // for normalised RDCs (Q1 = Rob Best's method, Q2 = Cornilescu method)
      cout << "# Q-values for back-calculated data" << endl;
      double Q1_norm = RDCTools.calc_Q1(calc_bc_norm, exp_bc_norm);
      fprintf(stdout, "# Normalized Q (R.Best) = %12.6f\n",Q1_norm);
      double Q2_norm = RDCTools.calc_Q2(calc_bc_norm, exp_bc_norm);
      fprintf(stdout, "# Normalized Q (PALES)  = %12.6f\n", Q2_norm);
      // for un-normalised RDCs

      double Q1_raw = RDCTools.calc_Q1(calc_bc, exp_bc);
      fprintf(stdout, "# Raw Q (R.Best) = %12.6f\n", Q1_raw);
      double Q2_raw = RDCTools.calc_Q2(calc_bc, exp_bc);
      fprintf(stdout, "# Raw Q (PALES)  = %12.6f\n", Q2_raw);
      // R-values
      cout << "# R-values for back-calculated data" << endl;
      double R_norm = RDCTools.calc_R(calc_bc_norm, exp_bc_norm);
      fprintf(stdout, "# Normalized R = %12.6f\n",R_norm);
      double R_raw = RDCTools.calc_R(calc_bc, exp_bc);
      fprintf(stdout, "# Raw R = %12.6f\n", R_raw);
      // rmsd
      // back-conversion factor to scale from ps-1 into Hz (s-1)
      double rdc_bcf = 1.0 / gmath::physConst.get_pico();
      cout << "# RMSD for back-calculated data" << endl;
      double RMSD_raw = RDCTools.calc_RMSD(calc_bc, exp_bc)*rdc_bcf;
      fprintf(stdout, "# Raw RMSD = %12.6f Hz\n", RMSD_raw);

      // -- print the experimental and back-calculated RDCs --
      cout.setf(ios::fixed, ios::floatfield);
      cout.precision(5);
      // title
      cout << setw(18) << "# residue/atom 1  " << setw(18) << "    residue/atom 2"
      << setw(18) << "    residue/atom 3" << setw(11) << "D_exp" << setw(11) << "D_calc" << endl;

      // now print actual (back-calculated) RDCs
      for (unsigned int i = 0; i < nbc; i++) {

        double calc_rdc = gsl_vector_get(calc_bc, i) * rdc_bcf;
        // k may be zero if not a side-chain NH RDC
        int katom, kres;
        string kname;
        if (bc_dat[i].type == 5 || bc_dat[i].type == 6) {
          kres = sys.mol(bc_dat[i].mol).topology().resNum(bc_dat[i].k) + 1;
          katom = bc_dat[i].k + 1;
          kname = sys.mol(bc_dat[i].mol).topology().atom(bc_dat[i].k).name();
        } else {
          katom = 0;
          kname = "--";
          kres = 0;
          // abs value for HH RDCs
          if (bc_dat[i].type == 7 || bc_dat[i].type == 8) {
            calc_rdc = std::abs(calc_rdc);
          }
        }

        cout << setw(6) << sys.mol(bc_dat[i].mol).topology().resNum(bc_dat[i].i) + 1
                << setw(6) << sys.mol(bc_dat[i].mol).topology().atom(bc_dat[i].i).name()
                << setw(6) << bc_dat[i].i + 1
                << setw(6) << sys.mol(bc_dat[i].mol).topology().resNum(bc_dat[i].j) + 1
                << setw(6) << sys.mol(bc_dat[i].mol).topology().atom(bc_dat[i].j).name()
                << setw(6) << bc_dat[i].j + 1
                << setw(6) << kres
                << setw(6) << kname
                << setw(6) << katom
                << setw(11) << gsl_vector_get(exp_bc, i) * rdc_bcf
                << setw(11) << calc_rdc
                << endl;
      }

      // clean up
      gsl_matrix_free(coef_mat);
      gsl_vector_free(exp_fit_norm);
      gsl_vector_free(exp_bc);
      gsl_vector_free(calc_bc);
      gsl_vector_free(calc_bc_norm);      
      if (backcalc) {
        gsl_vector_free(exp_bc_norm);
        gsl_matrix_free(coef_mat_bc_j); 
        gsl_matrix_free(coef_mat_bc_k); 
      }
    }

    // start JVAL-specific loop
    else if (datatype == "JVAL") {

      // parse boundary conditions for sys
      Boundary *pbc = BoundaryParser::boundary(sys, args);
      // parse gather method
      Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

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

      // declare all variables
      int nfit, nbc;
      vector<JvalData::jvalues> fit_dat;
      vector<JvalData::jvalues> bc_dat;

      // note Jval are not scaled (cf RDCs)
      gsl_vector *exp_fit, *exp_bc, *calc_bc, *K, *Ktmp, *work;
      gsl_matrix *coef_mat, *coef_mat_bc, *A, *V;

      // read in the Jvalue data
      PropertyContainer fit_props(sys, pbc);
      double delta;
      JvalTools.read_jval(buffer, sys, fit_dat, fit_props, delta);
      nfit = fit_dat.size();
      sf.close();

      // allocate sizes
      coef_mat = gsl_matrix_alloc(nfit, 3);
      gsl_matrix_set_zero(coef_mat);

      A = gsl_matrix_alloc(nfit, 3);
      V = gsl_matrix_alloc(3, 3);
      K = gsl_vector_alloc(3);
      Ktmp = gsl_vector_alloc(3);
      work = gsl_vector_alloc(3);
      exp_fit = gsl_vector_alloc(nfit);

      // now we can fill jval_vec
      JvalTools.fill_jvalvec(fit_dat, exp_fit, nfit);

      // check if we have different data to back-calculate
      PropertyContainer bc_props(sys, pbc);
      bool backcalc = false;
      if (args.count("calcspec") > 0) {
        backcalc = true;
        Ginstream sf(args["calcspec"]);
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
     
        JvalTools.read_jval(buffer, sys, bc_dat, bc_props, delta);
        nbc = bc_dat.size();

        sf.close();
      }// otherwise copy fit data into bc data
    else {
      bc_dat = fit_dat;
      nbc = bc_dat.size();
      bc_props = fit_props;
    }

    // initialise arrays for back-calculated JVAL data
    calc_bc = gsl_vector_alloc(nbc);
    exp_bc = gsl_vector_alloc(nbc);
    coef_mat_bc = gsl_matrix_calloc(nbc, 3);

    if (backcalc) {
      // fill exp_bc_sc with experimental Jval that will be back-calculated
      JvalTools.fill_jvalvec(bc_dat, exp_bc, nbc);
    }

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

          // and for the RDCs to back-calculate (if different)
          if (backcalc) {
            JvalTools.calc_coef(sys, bc_props, coef_mat_bc, nbc, w, delta);
          }

        }//end if use_this_frame
        if (done)
          break;

      }//end frame
    }// end io

    // scale coefficients by number of frames (sum of weights as each contribution
    // was multiplied by its weight)
    gsl_matrix_scale(coef_mat, 1. / weight_sum);

    // if we are back-calculating the same data that we fit to, we can just
    // copy the coefficients and RDCs across
    if (!backcalc) {
      // back-calculating same data as fitted to so just copy coefficients
      coef_mat_bc = coef_mat;
      exp_bc = exp_fit;
    }// if not, the coefficients were already calculated and exp_bc already filled
        // so all we have to do is scale the coefficients
      else {
        gsl_matrix_scale(coef_mat_bc, 1. / weight_sum);
      }

      // make a copy of coef_mat because it is modified by the SVD functions
      gsl_matrix_memcpy(A, coef_mat);

      // calculate and print Karplus parameters (fitting to fit data)
      gsl_linalg_SV_decomp(A, V, Ktmp, work);
      gsl_linalg_SV_solve(A, V, Ktmp, exp_fit, K);

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
      fprintf(stdout, "# A = %12.5e\n", K_A);
      fprintf(stdout, "# B = %12.5e\n", K_B);
      fprintf(stdout, "# C = %12.5e\n", K_C);

      // back-calculate J-values
      gsl_blas_dgemv(CblasNoTrans, 1.0, coef_mat_bc, K, 0., calc_bc);
      // calculate and print Q
      double Q = JvalTools.calc_Q(calc_bc, exp_bc);
      cout << "# Goodness of fit:" << endl;
      fprintf(stdout, "# Q = %12.6f\n", Q);

      cout.setf(ios::fixed, ios::floatfield);
      cout.precision(5);

      // title for back-calculated Jvalues
      cout << setw(5) << "#   i" << setw(5) << "j" << setw(5) << "k" << setw(5) << "l"
              << setw(11) << "J_calc" << setw(11) << "J_exp" << endl;

      // now print actual (back-calculated) Jvalues
      for (int i = 0; i < nbc; i++) {

        cout << setw(5) << bc_dat[i].i + 1
                << setw(5) << bc_dat[i].j + 1
                << setw(5) << bc_dat[i].k + 1
                << setw(5) << bc_dat[i].l + 1
                //<< setw(11) << gsl_vector_get(bc, i) * bcf
                //<< setw(11) << gsl_vector_get(jval_vec_bc, i) * bcf
                << setw(11) << gsl_vector_get(calc_bc, i)
                << setw(11) << gsl_vector_get(exp_bc, i)
                << endl;
      }

      // clean up??
      gsl_vector_free(exp_fit);
      gsl_matrix_free(coef_mat);

    } else {
      throw gromos::Exception("svd_fit", "Data types other than RDC or JVAL"
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

