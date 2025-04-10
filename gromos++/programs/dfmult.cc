/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file dfmult.cc
 * multiple free energy differences
 */

/**
 * @page programs Program Documentation
 *
 * @anchor dfmult
 * @section dfmult multiple free energy differences 
 * @author @ref cc
 * @date 24. 04. 09
 *
 * Calculates free energy differences between multiple states @f$A@f$ and
 * @f$B@f$ from a simulation at a reference state @f$R@f$ according to
 * @f[ \Delta F_{BA} = F_B - F_A
 *                   = F_B - F_R - (F_A - F_R)
 *                   = \Delta F_{BR} - \Delta F_{AR}
 *                   = -\beta^{-1} \ln \frac{\left<\exp \left[-\beta \left( V_B -V_R\right) \right]\right>_R}{\left<\exp \left[-\beta \left( V_A -V_R\right) \right]\right>_R}@f]
 * It reads in energy time series generated by @ref ene_ana. It provides an  error
 * estimate (err) which is based on calculation of the
 * (co)variances and the statistical inefficiency as described in
 * J. Chem. Theory Comput. 2007, 3, 26-41.
 * The implementation closely follows the Python implementation provided by the
 * authors. When calculating averages and uncertainties special care is taken
 * in order to avoid overflow (see Comput. Phys. Comm. 2003, 153, 397-406).
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@temp</td><td>&lt;temperature for perturbation&gt; </td></tr>
 * <tr><td> \@stateR</td><td>&lt;energy file for state R&gt; </td></tr>
 * <tr><td> \@endstates</td><td>&lt;energy files of end states&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
 dfmult
    @temp      300
    @stateR    eR.dat
    @endstates e1.dat e2.dat e3.dat e4.dat
    @endverbatim
 *
 * <hr>
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Physics.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;


int main(int argc, char** argv) {
  Argument_List knowns;

  knowns << "temp" << "stateR" << "endstates";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@temp      <temperature for perturbation>\n";
  usage +=   "\t@stateR    <energy file for state R>\n";
  usage +=   "\t@endstates <energy files of endstates>\n";
  
  try {
    Arguments args(argc, argv, knowns, usage);

    // Get temperature as a double
    double temp = args.getValue<double>("temp");

    // store all the time series for the error analysis (BIG!!)
    vector<gmath::Stat< double> > allexphxhr;

    // make sure we have the @stateR flag
    args.check("stateR",1);
    
    // loop over endstates
    args.check("endstates", 1);
    for (Arguments::const_iterator
      iter = args.lower_bound("endstates"),
            to = args.upper_bound("endstates");
            iter != to; ++iter) {

      // open end state energy file
      ifstream stateX;
      stateX.open((iter->second).c_str());
      if (!stateX)
        throw gromos::Exception("dfmult", "Could not open energy file for end state\n");

      // open reference state energy file
      ifstream stateR;
      Arguments::const_iterator iter2 = args.lower_bound("stateR");

      stateR.open((iter2->second).c_str());
      if (!stateR)
        throw gromos::Exception("dfmult", "Could not open energy file for state R\n");

      // check for end of file
      bool eofR = false, eofX = false;
      // check for reading errors
      bool errorR = false, errorX = false;
      // read in the time (later perhaps add @tmin and @tmax flags)
      double timeR, timeX;
      // variables to store the energies;
      double hr, hx;
      
      string sdum;
      // save the -b(H_X - H_R) values
      gmath::Stat< double> exphxhr;
      
      while (true) {
        // read from state R
        while (true) {
          std::getline(stateR, sdum);
          if (stateR.eof()) {
            eofR = true;
            break;
          }
          std::string::size_type it = sdum.find('#');
          if (it != std::string::npos)
            sdum = sdum.substr(0, it);
          if (sdum.find_first_not_of(" \t") == std::string::npos)
            continue;
          std::istringstream is(sdum);
          if (!(is >> timeR >> hr))
            errorR = true;
          break;
        }

        // read from endstate (this could be taken out of the loop over endstates)
        while (true) {
          std::getline(stateX, sdum);
          if (stateX.eof()) {
            eofX = true;
            break;
          }
          std::string::size_type it = sdum.find('#');
          if (it != std::string::npos)
            sdum = sdum.substr(0, it);
          if (sdum.find_first_not_of(" \t") == std::string::npos)
            continue;
          std::istringstream ix(sdum);
          if (!(ix >> timeX >> hx))
            errorX = true;
          break;
        }

        // check eof / error
        if (eofR && eofX) break;
        if (eofX || errorX || eofR || errorR)
          throw gromos::Exception("dfmult", "Could not read from energy file. Check lengths!");

        exphxhr.addval((-(hx - hr) / (gmath::physConst.get_boltzmann_silent() * temp)));  
      }

      stateX.close();
      stateR.close();

      // store the time series for error analysis
      allexphxhr.push_back(exphxhr);
      
    } // end loop over end states

    //cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::scientific, ios::floatfield);
    cout.precision(7);
    cout << setw(36) << "#DF (kJ/mol)" << setw(18) << "err"  << endl;
    // loop over all saved values and calculate free energy differences and uncertainties
    for(unsigned int i=0;i<allexphxhr.size();i++) {
      //cout << "# i: " << i + 1 << endl;
      // number of data points
      int N = allexphxhr.at(i).n();
      // df
      double df_ir = allexphxhr.at(i).lnexpave();
      // variance
      int sign_var_ii = 1;
      double var_ii = gmath::Stat<double>::lnexpcovariance(allexphxhr.at(i), allexphxhr.at(i), sign_var_ii);
      if(sign_var_ii < 0)
        throw gromos::Exception("dfmult", "Got negative variance!");
      // statistical inefficiency
      double si_ii = gmath::Stat<double>::lnexp_stat_ineff(allexphxhr.at(i), allexphxhr.at(i));
      // error estimate contribution
      double d2i = var_ii - log(double(N)) + si_ii;
      stringstream name;
      name << "DF_" << i + 1 << "_R";
      cout << setw(18) << name.str().c_str()
              << setw(18) << -gmath::physConst.get_boltzmann_silent() * temp * df_ir
              << setw(18) << gmath::physConst.get_boltzmann_silent() * temp * sqrt(exp(d2i - 2 * df_ir)) << endl;

      for(unsigned int j=(i+1);j<allexphxhr.size();j++){

        //cout << " # j: " << j + 1 << endl;
        // df
        double df_jr = allexphxhr.at(j).lnexpave();
        // variance
        int sign_var_jj = 1;
        double var_jj = gmath::Stat<double>::lnexpcovariance(allexphxhr.at(j), allexphxhr.at(j), sign_var_jj);
        if(sign_var_jj < 0)
          throw gromos::Exception("dfmult","Got negative variance!");
        // statistical inefficiency
        double si_jj = gmath::Stat<double>::lnexp_stat_ineff(allexphxhr.at(j), allexphxhr.at(j));
        // error estimate contribution
        double d2j = var_jj - log(double(N)) + si_jj;

        // df
        double df_ji = allexphxhr.at(j).lnexpave() - allexphxhr.at(i).lnexpave();
        // variance
        int sign_var_ji = 1;
        double var_ji = gmath::Stat<double>::lnexpcovariance(allexphxhr.at(j), allexphxhr.at(i), sign_var_ji);
        // statistical inefficiency
        double si_ji = gmath::Stat<double>::lnexp_stat_ineff(allexphxhr.at(j), allexphxhr.at(i));
        // Ensure that estimate of cross-uncertainty is reasonable (see Appendix).
        double max_si_ji = log(sqrt(exp(si_ii + si_jj))
                * abs(sqrt(exp(var_ii + var_jj - 2 * var_ji))) );
        // log is strictly monotone so < and > holds also for the logarithms
        if ((si_ji > max_si_ji) && (max_si_ji > log(1.0))) {
          cerr << "setting si_ji to max_si_ji" << endl;
          si_ji = max_si_ji;
        }
        
        // error estimate contribution
        double djdi = var_ji - log(double(N)) + si_ji;
        stringstream name;
        name << "DF_" << j + 1 << "_" << i + 1;
        cout << setw(18) << name.str().c_str()
                << setw(18) << -gmath::physConst.get_boltzmann_silent() * temp * df_ji
                << setw(18) << gmath::physConst.get_boltzmann_silent() * temp *
                               sqrt(exp(d2i - 2 * df_ir)
                                  + exp(d2j - 2 * df_jr)
                              - 2 * exp(djdi - (df_jr + df_ir)) * sign_var_ji)
                << endl;
      } // states j
    } // states i


  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

