// utils_RdcFuncs.cc

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "RdcFuncs.h"
#include "../args/Arguments.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Physics.h"
#include "../gmath/Vec.h"
#include "../gromos/Exception.h"
#include "../gio/Ginstream.h"

using namespace std;
using namespace gio;
using namespace gcore;
using namespace gmath;

using gcore::System;
using args::Arguments;
using utils::RdcFuncs;

// Constructor
RdcFuncs::RdcFuncs(System &sys, Arguments &args) {}

// Destructor
RdcFuncs::~RdcFuncs() {}

// function to read RDC data
void RdcFuncs::read_rdc(vector<string> buffer, const System &sys, vector<RDCData::rdcparam> &rdcp) {

  // read into rdcparam class
  for (unsigned int jj = 1; jj < buffer.size() - 1; jj++) {

    // first get atom numbers
    istringstream is(buffer[jj]);
    int i, j;
    is >> i >> j;

    // adjust i and j to gromos numbering
    i--;
    j--;

    // define local set of RDC parameters   
    RDCData::rdcparam rp;

    // get molecule number, offset atoms numbers
    int m = 0, offset = 0;
    for (; i >= sys.mol(m).numAtoms() + offset; m++)
      offset += sys.mol(m).numAtoms();
    rp.mol = m;
    rp.i = i - offset;
    rp.j = j - offset;

    // now get the other parameters
    is >> rp.w >> rp.exp >> rp.gi >> rp.gj;

    // convert into gromos units
    double rdc_cf = pico;
    double g_cf = ( atomic_mass_unit / elementary_charge ) * 1.0e7;
    rp.exp *= rdc_cf;
    rp.gi *= g_cf;
    rp.gj *= g_cf;

    if (is.fail())
      throw gromos::Exception("svd_fit", "Bad line in rdc-file\n" + buffer[jj]);

    rdcp.push_back(rp);
  }

}

// function to compute rij and Dmax from first frame
void RdcFuncs::calc_dmax(const System &sys, vector<RDCData::rdcparam> &rdcp,
        int nrdc) {

  for (int d = 0; d < nrdc; d++) {

    // get the distance rij
    double rij = (sys.mol(rdcp[d].mol).pos(rdcp[d].i) - sys.mol(rdcp[d].mol).pos(rdcp[d].j)).abs();

    // compute and store Dmax
    rdcp[d].dmax = (mu0 * h * rdcp[d].gi * rdcp[d].gj) / (8.0 * pi * pi * pi * rij * rij * rij);
    
  }
}

// function to read weights for individual frames from file
void RdcFuncs::read_weights(vector<string> buffer, vector<RDCWeights::weights> &weight_data) {

  // read into weights vector
  for (int d = 1; d < buffer.size() - 1; d++) {

    // define local weight struct
    RDCWeights::weights w;
    // get values from file
    istringstream is(buffer[d]);
    is >> w.frame >> w.weight;
    // add to weight_data
    weight_data.push_back(w);

  }
}

// compute the coefficients of the matrix describing bond vector fluctuations
void RdcFuncs::calc_coef(const System &sys, vector<RDCData::rdcparam> &fit_data,
        gsl_matrix *coef_mat, int nrdc, double w)
{
  int i, j, m;
  double mu_x, mu_y, mu_z, mu_r;

  for (int d = 0; d < nrdc; d++) {

    // get atom numbers
    i = fit_data[d].i;
    j = fit_data[d].j;
    m = fit_data[d].mol;

    // get diff in coords
    mu_x = sys.mol(m).pos(i)[0] - sys.mol(m).pos(j)[0];
    mu_y = sys.mol(m).pos(i)[1] - sys.mol(m).pos(j)[1];
    mu_z = sys.mol(m).pos(i)[2] - sys.mol(m).pos(j)[2];

    mu_r = sqrt(mu_x * mu_x + mu_y * mu_y + mu_z * mu_z);
    mu_x /= mu_r;
    mu_y /= mu_r;
    mu_z /= mu_r;

    gsl_matrix_set(coef_mat, d, 0,
            gsl_matrix_get(coef_mat, d, 0) + w * (mu_x * mu_x - mu_z * mu_z));
    gsl_matrix_set(coef_mat, d, 1,
            gsl_matrix_get(coef_mat, d, 1) + w * (mu_y * mu_y - mu_z * mu_z));
    gsl_matrix_set(coef_mat, d, 2,
            gsl_matrix_get(coef_mat, d, 2) + w * (2.0 * mu_x * mu_y));
    gsl_matrix_set(coef_mat, d, 3,
            gsl_matrix_get(coef_mat, d, 3) + w * (2.0 * mu_x * mu_z));
    gsl_matrix_set(coef_mat, d, 4,
            gsl_matrix_get(coef_mat, d, 4) + w * (2.0 * mu_y * mu_z));
  }
}

// fill a vector with the normalised RDCs (scaled by Dmax)
void RdcFuncs::fill_rdcvec(const vector<RDCData::rdcparam> &R, gsl_vector *v, int nrdc)
{
	for (int i=0; i<nrdc; i++) {
		gsl_vector_set(v,i,R[i].exp/R[i].dmax);
	}
}

// calculate the Q value (goodness of fit)
double RdcFuncs::calc_Q(gsl_vector *calc, gsl_vector *expt)
{
	int ndat;
	double sdev2,scalc2,tmpf,expt_i,calc_i;
	sdev2 = 0;
	scalc2 = 0;
	ndat = calc->size;
	for (int i=0; i<ndat; i++) {
		calc_i = gsl_vector_get(calc,i);
		expt_i = gsl_vector_get(expt,i);
		tmpf = calc_i-expt_i;
		sdev2 += tmpf*tmpf;
		scalc2 += calc_i*calc_i;
	}
	return sqrt(sdev2/scalc2);
}

// fill a vector with un-normalised RDCs
void RdcFuncs::unnorm_rdcvec(const vector<RDCData::rdcparam> &R, gsl_vector *v, gsl_vector *w)
{
	int nrdc = R.size();
	double tmp;
	for (int i=0; i<nrdc; i++) {
		tmp = gsl_vector_get(v,i);
		gsl_vector_set(w,i,tmp*R[i].dmax);
	}
}
