/**
 * @file rdc_sh.cc
 * find probability distribution of magnetic field vector directions compatible
 * with a set of experimental RDCs by fitting to spherical harmonics
 */
/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor rdc_sh
 * @section rdc_sh find probability distribution in terms of spherical harmonics
 * @author @ref ja
 * @date 05.03.2010
 *
 * PROGRAM DESCRIPTION
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@rdcspec</td><td>&lt;file containing RDC to fit to&gt; </td></tr>
 * <tr><td> \@lsh</td><td>&lt;l for spherical harmonics&gt; </td></tr>
 * <tr><td> \@coord</td><td>&lt;coordinate file&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   rdc_sh
     @topo       ex.top
     @pbc        r
     @rdcspec    ex.rdc
     @lsh        3
     @coord      ex.g96
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
#include <complex>
#include <cmath>


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_legendre.h>
#include<gsl/gsl_multifit.h>

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
#include "../src/utils/PropertyContainer.h"


using namespace args;
using namespace bound;
using namespace fit;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace std;
using namespace utils;

// factorial function
int factorial (int num);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "rdcspec" << "lsh" << "coord";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>]\n";
  usage += "\t@rdcspec     <file containing RDC data to fit to>\n";
  usage += "\t@lsh         <l parameter for spherical harmonics>\n";
  usage += "\t@coord       <file containing coordinates>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // first print a warning
    cout << "# WARNING: this program is experimental. Please consult Jane before using\n" << endl;

    // define physical constants
    const double degree2radian = gmath::physConst.get_degree2radian();
    const double pi = gmath::physConst.get_pi();
    const double pico = gmath::physConst.get_pico();

    // read topology and initialise systems
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // define input coordinate
    InG96 ic;
    // open coordinate file
    ic.open(args["coord"]);
    ic.select("SOLUTE");
    // read coordinates
    ic >> sys;
    if (!sys.hasPos)
      throw gromos::Exception("rdc_sh",
            "Unable to read POSITION(RED) block from "
            "coordinate file.");
    // gather
    (*pbc.*gathmethod)();

    // access RDC-specific io and fitting functions
    RdcFuncs RDCTools(sys, args);

    // get experimental RDCs from file (in Hz) (also gyr)
    Ginstream sf(args["rdcspec"]);
    vector<string> buffer;
    vector<RDCData::rdcparam> rdc_exp;
    // search for RDCVALRESSPEC block
    bool found_rdc = false;
    while (!sf.stream().eof() && !found_rdc) {
      sf.getblock(buffer);
      if (buffer[0] == "RDCVALRESSPEC")
        found_rdc = true;
    }
    if (!found_rdc)
      throw gromos::Exception("rdc_sh", "RDC file does not contain a RDCVALRESSPEC block!");
    // check for END
    if (buffer[buffer.size() - 1].find("END") != 0)
      throw gromos::Exception("rdc_sh", "RDC file is corrupted. No END in "
            + buffer[0] + " block. Got\n" + buffer[buffer.size() - 1]);
    // read in the RDC data
    string calc_rij = "SPEC"; // HACK!!
    RDCTools.read_rdc(buffer, sys, rdc_exp, calc_rij, true);
    unsigned int nrdc = rdc_exp.size();
    sf.close();

    // read lsh
    vector<int> ltmp;
    if (args.count("lsh") > 0) {

      for (Arguments::const_iterator it = args.lower_bound("lsh");
              it != args.upper_bound("lsh"); ++it) {
        int bla = atoi(it->second.c_str());
        ltmp.push_back(bla);
      }
      if (ltmp.size() != 1) {
        throw gromos::Exception("rdc_sh",
                "you can only give one value for @lsh");
      }
    } else {
      throw gromos::Exception("rdc_sh",
              "you must specify a value for lsh");
    }
    const int lsh = ltmp[0];
    // compute total number of spherical harmonic functions
    unsigned int nsh = 0;
    for (int l = 0; l <= lsh; ++l) {
      nsh += 2 * l + 1;
    }
    if (lsh > int(nrdc)) {
      throw gromos::Exception("rdc_sh",
              "lsh cannot be larger than the number of RDCs");
    }

    // gsl declarations
    gsl_vector *rdc_exp_sc, *rdc_tmp, *Ylm_coeff;
    gsl_matrix *rdc_Ylm, *cov;

    // gsl allocations
    // scaled experimental rdcs
    rdc_exp_sc = gsl_vector_alloc(nrdc);
    // temporary storage for calculated rdcs
    rdc_tmp = gsl_vector_alloc(nrdc);
    // Ylm for +ve or -ve m
    //Ylm_posm = gsl_vector_alloc(nsh);
    //Ylm_negm = gsl_vector_alloc(nsh);
    // RDC * Ylm (zeroed)
    rdc_Ylm = gsl_matrix_calloc(nrdc, nsh);


    // THIS IS SOMEWHAT POINTLESS?? JUST HAVE TO PUT IT BACK LATER...
    // compute dmax (i.e. prefactor)
    RDCTools.calc_dmax16(rdc_exp);
    // put scaled rdcs (/dmax) into gsl vector
    //RDCTools.fill_rdcvec_norm(rdc_exp, rdc_exp_sc);

    // DEBUG: try without scaling them
    for (unsigned int i = 0; i < nrdc; i++) {
      gsl_vector_set(rdc_exp_sc, i, rdc_exp[i].exp);
    }

    // DEBUG
    //cout << rdc_exp[0].i << " " << rdc_exp[0].j << " RDCexp: " << rdc_exp[0].exp << " gi: "
    //<< rdc_exp[0].gi << " gj: " << rdc_exp[0].gi << " Dmax: " << rdc_exp[0].dmax <<
    //        " scaledRDC: " << gsl_vector_get(rdc_exp_vec,0) << endl;


    // loop over theta and phi (in deg for now)
    for (unsigned int theta = 0; theta < 90; ++theta) {

      // put theta into radians and compute cos(theta)
      double theta_rad = theta * degree2radian;
      double cos_theta_rad = cos(theta_rad);

      // uniform distribution
      const double uniform = (1. / (2. * pi) ) * sin(theta_rad);

      // loop over phi (in deg for now)
      for (unsigned int phi = 0; phi < 360; ++phi) {

        // DEBUG
        //cout << "\n-----------theta = " << theta << "-------------\n" << endl;
        //cout << "\n-------------phi = " << phi << "-------------\n" << endl;

        // initiate counter of spherical harmonic functions
        unsigned int this_sh = 0;

        // convert phi to radians
        double phi_rad = phi * degree2radian;

        // compute RDCs for this orientation of H
        RDCTools.calc_rdc_H(sys, rdc_exp, theta_rad, phi_rad, rdc_tmp);

        // loop over values of l
        for (int l = 0; l <= lsh; ++l) {

          // first do m = 0
          // get normalised Legendre polynomials P^l_m
          double Pl_m0 = gsl_sf_legendre_sphPlm(l, 0, cos_theta_rad);
          //double Pl_m0 = gsl_sf_legendre_Plm(l,0,cos_theta_rad);
          // multiply P^l_m by exp(i*m*phi) to get Y^l_m
          // exp(i*m*phi) always 1 if m = 0; this is pointless, just for naming compatibility
          double Yl_m0 = Pl_m0;

          // multiply each RDC by Y^l_m and store (sum) in rdc_Ylm matrix
          RDCTools.rdcYlm(nrdc, rdc_tmp, this_sh, Yl_m0, rdc_Ylm);


          // and now for m > 0 (m must be <= l)
          for (int m = 1; m <= l; ++m) {

            // first increment this_sh
            ++this_sh;

            // compute normalised P^l_m for +ve m
            double Pl_posm = gsl_sf_legendre_sphPlm(l, m, cos_theta_rad);
            //double Pl_posm = gsl_sf_legendre_Plm(l,m,cos_theta_rad);
            // scale to get P^l_m for -ve m
            double denominator = factorial(l+m); // doesn't like two factorial calls below
            double Pl_negm = pow(-1.0,double(m)) * ( (factorial(l-m)) / denominator) * Pl_posm;
            //double Pl_negm = pow(-1.0, double(m)) * ((factorial(l - m)) / (factorial(l + m))) * Pl_posm;

            // DEBUG
            //cout << "l = " << l << " m = " << m << " -1^m = " << pow(-1,m) <<
            //" (l-m)! = " << factorial(l-m) << " (l+m)! = " << factorial(l+m) << endl;
            //cout << "Pl+m = " << Pl_posm << "\tPl-m = " << Pl_negm << endl;

            
            // compute exp(i*m*phi)
            complex<double> pow_iposmphi(0., double(m) * phi_rad);
            complex<double> pow_inegmphi(0., double(-m) * phi_rad);
            complex<double> exp_iposmphi = exp(pow_iposmphi);
            complex<double> exp_inegmphi = exp(pow_inegmphi);
            // multiply P^l_m by exp(i*m*phi) to get Y^l_m
            double Yl_posm = Pl_posm * exp_iposmphi.imag();
            double Yl_negm = Pl_negm * exp_inegmphi.imag();
             

            //double Yl_posm = Pl_posm * cos(m * phi_rad);
            //double Yl_negm = Pl_negm * cos(m * phi_rad);

            // DEBUG
            //cout << "Complex power +m real = " << pow_iposmphi.real() <<
            //" and imag = " << pow_iposmphi.imag() << "\nComplex power -m real = "
            //<< pow_inegmphi.real() << " and imag = " << pow_inegmphi.imag() << endl;
            //cout << "e(i+mphi).real = " << exp_iposmphi.real() << " and imag = "
            //<< exp_iposmphi.imag() << "\ne(i-mphi).real = " << exp_inegmphi.real()
            //<< " and imag = " << exp_inegmphi.imag() << endl;

            // multiply each RDC by Y^l_m and store (sum) in rdc_Ylm matrix
            RDCTools.rdcYlm(nrdc, rdc_tmp, this_sh, Yl_posm, rdc_Ylm);

            // increment this_sh before dealing with -ve m
            ++this_sh;
            RDCTools.rdcYlm(nrdc, rdc_tmp, this_sh, Yl_negm, rdc_Ylm);

          } // end m loop

          // increase spherical harmonic counter
          ++this_sh;

        } // end l loop
      } // end phi loop

      // multiply stored rdc*Ylm by uniform distribution for this theta
      gsl_matrix_scale(rdc_Ylm, uniform);

    } // end theta loop

    // DEBUG: print some values of rdc_Ylm for l=1,m=1
    //for (int i = 0; i < nrdc; ++i) {
    //  cout << "rdc " << i+1 << "\trdc_Ylm " << gsl_matrix_get(rdc_Ylm, i, 2) << endl;
    //}

    // CURRENTLY:
    // rdc_Ylm contains sum("integral") of RDCcalc (without prefactor) * Y^l_m
    // for each RDC (internuclear vector)[n] and for each combination of l and m [p]
    // rdc_Ylm[n,p]*COEFF[p] = rdc_exp_sc[n]: least-squares to estimate COEFF

    // definitions and allocations
    // workspace for fitting
    gsl_multifit_linear_workspace * work;
    work = gsl_multifit_linear_alloc(nrdc, nsh);
    // chi-square returned by gsl least-squares fit
    double chisq;
    // covariance matrix returned by gsl least-squares fit
    cov = gsl_matrix_alloc(nsh, nsh);
    // fitted coefficients of the spherical harmonics
    Ylm_coeff = gsl_vector_alloc(nsh);

    // call the least-squares fit
    gsl_multifit_linear(rdc_Ylm, rdc_exp_sc, Ylm_coeff, cov, &chisq, work);


    // the "svd" options allows you to specify the tolerance before singular
    // values are discarded and returns the rank (so you can see if any were discarded)
    //size_t rank;
    //gsl_multifit_linear_svd(rdc_Ylm, rdc_exp_sc, 1e-20, &rank, Ylm_coeff,
    //        cov, &chisq, work);

    /*
    // DEBUG: check derivative
    cout << "# Derivatives of least-squares fitting function:" << endl;
    double sum_dv_dc = 0;
    for (int l = 0; l < nsh; ++l) {
      double dv_dc = 0;
      for (int r = 0; r < nrdc; ++r) {
        for (int ll = 0; ll < nsh; ++ll) {
          dv_dc += (gsl_matrix_get(rdc_Ylm, r, ll) * gsl_vector_get(Ylm_coeff, ll));
        }
        dv_dc -= gsl_vector_get(rdc_exp_sc, r); // 1e-4
        dv_dc *= gsl_matrix_get(rdc_Ylm, r, l); // +/- 20
      }
      cout << "# sh = " << l << ", dV/dc_lm " << dv_dc << endl;
      sum_dv_dc += dv_dc;
    }
    cout << "# Sum of dV/dc_lm " << sum_dv_dc << endl;
     */

    // OUTPUT
    cout.precision(10);
    cout.setf(ios::right, ios::adjustfield);
    cout << "# Spherical harmonic coefficients" << endl;
    cout << "#" << setw(5) << "l" << setw(5) << "m" << setw(14) << "coeff" << endl;
    unsigned int mm = 0;
    for (int l = 0; l <= lsh; ++l) {
      // for m = 0
      fprintf(stdout, "# %4d %4d %13.5e\n", l, 0, gsl_vector_get(Ylm_coeff, mm));

      for (int m = 1; m <= l; ++m) {

        // for +m
        ++mm;
        fprintf(stdout, "# %4d %4d %13.5e\n", l, m, gsl_vector_get(Ylm_coeff, mm));
        // for -m
        ++mm;
        fprintf(stdout, "# %4d %4d %13.5e\n", l, -m, gsl_vector_get(Ylm_coeff, mm));

      }
    }
    // chi square
    fprintf(stdout, "#\n# ChiSquare = %13.10f\n", chisq);
    // if svd, print rank
    //fprintf(stdout, "# Rank = %12d\n", rank);

    // title for actual RDCs
    cout << setw(9) << "#\n#   atom1" << setw(9) << "atom2"
            << setw(11) << "RDC_exp" << setw(18) << "RDC_calc" << endl;

    // compare rdc_Ylm * Tlm_coeff to experimental RDCs
    for (unsigned int d = 0; d < nrdc; d++) {
      double rdc_calc_tmp = 0;
      for (unsigned int mm = 0; mm < nsh; mm++) {

        // get contribution from each spherical harmonic to this RDC
        rdc_calc_tmp += gsl_matrix_get(rdc_Ylm, d, mm) * gsl_vector_get(Ylm_coeff, mm);
      }

      cout << setw(5) << rdc_exp[d].i + 1
              << setw(4) << sys.mol(rdc_exp[d].mol).topology().atom(rdc_exp[d].i).name()
              << setw(5) << rdc_exp[d].j + 1
              << setw(4) << sys.mol(rdc_exp[d].mol).topology().atom(rdc_exp[d].j).name()
              << setw(11) << rdc_exp[d].exp / pico
              // DEBUG: not scaled
              //<< setw(18) << (rdc_calc_tmp * rdc_exp[d].dmax) / pico
              << setw(18) << rdc_calc_tmp / pico
              << endl;

    }

    // free workspace
    gsl_multifit_linear_free(work);


  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

// factorial function
int factorial (int num)
{
 if (num<=1)
  return 1;
 return factorial(num-1)*num; // recursive call
}
