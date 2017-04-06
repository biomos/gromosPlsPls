/**
 * @file extendedTI.cc
 * reads energy and free energy trajectories and reweights
 * to specific properties
 */

/**
 * @page programs Program Documentation
 *
 * @anchor extendedTI
 * @section extendedTI reweight (free) energy trajectories
 * @author @ref adr
 * @date 18. 8. 2015
 *
 * INSERT DESCRIPTION HERE!
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@en_files</td><td>&lt;energy files&gt; (and/or) </td></tr>
 * <tr><td> \@fr_files</td><td>&lt;free energy files&gt; </td></tr>
 * <tr><td> \@library</td><td>&lt;library for property names&gt; [print] </td></tr>
 * <tr><td> \@prop</td><td>&lt;@ref predicting [lambda or NLAM]&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  extendedTI
    @en_files   ex.tre
    @fr_files   ex.trg
    @library    ene_ana.lib
    @prop	lambda

   @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <map>
#include <cmath>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/EnergyTraj.h"
#include "../src/gmath/Expression.h"

using namespace std;
using namespace args;
using namespace gio;
using namespace gcore;

void set_standards(utils::EnergyTraj &e);
void read_library(string name, utils::EnergyTraj& e);

int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "en_files" << "fr_files" << "library" << "nrlambdas" << "minlam" << "maxlam" 
         << "slam" << "plam" << "NLAMs" << "NLAMp" << "temp" << "no_error" << "printX";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@en_files    <energy files> and\n";
  usage += "\t@fr_files    <free energy files>\n";
  usage += "\t@temp        <temperature>\n";
  usage += "\t@library     <library for block information> \n";
  usage += "\t@nrlambdas   <nr precalculated lambdas>\n";
  usage += "\t@minlam      <minimum precalculate lambda>\n";
  usage += "\t@maxlam      <maximum precalculated lambda>\n";
  usage += "\t@slam        <lambda value of simulation>\n";
  usage += "\t[@plam       <only define if predicting specific lambda value]>\n";
  usage += "\t@NLAMs       <NLAM value of simulation>\n";
  usage += "\t@NLAMp       <NLAM value to predict for >\n";
  usage += "\t[@no_error   <error estimates are not calculated if this flag is set >\n";
  usage += "\t[@printX     <write out values for plam to X.dat>\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    // check whether we have both energy and free energy files
    if(args.count("en_files")<=0 || args.count("fr_files")<=0)
      throw gromos::Exception("extendedTI", "en_files and fr_files need to be specified:\n"+usage);

    // do we predict for a single lambda, or all?
    bool single_plam;
    double plam;
    if(args.count("plam")<=0) 
      single_plam = false;
    else {
      single_plam = true;
      plam = args.getValue<double>("plam", true);
    }

    // do we calculate error estimates?
    bool no_error; 
    if(args.count("no_error")<0) no_error = false;
    else no_error = true; 

    // do we print X
    bool printX;
    double print_plam;
    if(args.count("printX")<0) printX = false;
    else {
      printX = true;
      print_plam = args.getValue<double>("printX",false);
    }

    double temp = args.getValue<double>("temp", true);
    int nrlambdas = args.getValue<int>("nrlambdas", true);
    double minlam = args.getValue<double>("minlam", true);
    double maxlam = args.getValue<double>("maxlam", true);
    double slam = args.getValue<double>("slam", true);
    double NLAMs = args.getValue<double>("NLAMs", true);
    double NLAMp = args.getValue<double>("NLAMp", true);

    
    // NEW: require a library 
    if(args.count("library") <=0)
      throw gromos::Exception("extendedTI", "no library file specified:\n"+usage);
    
    // read a library file?
    string library="";
    {
      Arguments::const_iterator iter=args.lower_bound("library"), 
	to=args.upper_bound("library");
      if(iter!=to){
	library=iter->second;
	++iter;
      }
    }

    // determine the index corresponding to slam 
    double lam_step = (maxlam - minlam)/(nrlambdas-1);
    int sindex = int(round(slam/lam_step));
    int pindex = int(round(plam/lam_step));
    int nr_plam; 

    if(single_plam){
      nr_plam = 1;
    }
    else {
      nr_plam = nrlambdas; 
    }

    // declare some vectors which we will need later
    vector<gmath::Stat<double> > x(nr_plam);
    vector<gmath::Stat<double> > vyvr(nr_plam);
    vector<gmath::Stat<double> > expvyvr(nr_plam);
    vector<gmath::Stat<double> > Xexpvyvr(nr_plam);

    // define an energy trajectory
    utils::EnergyTraj etrj;

    // read topology for the mass and the number of molecules
    double mass=0;
    etrj.addConstant("MASS",mass);
    
    // learn about the variable names how they map to the elements
    read_library(library, etrj);
    
    // define two input streams
    Ginstream gin_en;
    Ginstream gin_fr;
    bool do_energy_files     =(args.count("en_files")>0);
    bool do_free_energy_files=(args.count("fr_files")>0);
    
    Arguments::const_iterator it_en=args.lower_bound("en_files"),
      to_en=args.upper_bound("en_files"),
      it_fr=args.lower_bound("fr_files"),
      to_fr=args.upper_bound("fr_files");
    int cont=0, en_cont=0, fr_cont=0;
    if(do_energy_files) {
      gin_en.open(it_en->second.c_str()); 
      ++it_en; 
      en_cont=1;
    }
    
    if(do_free_energy_files) {
      gin_fr.open(it_fr->second.c_str());
      ++it_fr;
      fr_cont=1;
    }
    
    cont=en_cont+fr_cont;

    // starting with the calculations    
    while (true) {
      // read the next frame of the energy trajectory
      if(do_energy_files){
	int end_of_file=etrj.read_frame(gin_en, "ENERTRJ");
	if(end_of_file){
	  if(it_en!=to_en){
	    gin_en.close();
	    gin_en.open(it_en->second.c_str());
	    ++it_en;
	    //try again...
	    etrj.read_frame(gin_en, "ENERTRJ");
	  }
	  else
	    break;
	}
      }
      // read the next frame of the free energy trajectory
      if(do_free_energy_files){
	int end_of_file=etrj.read_frame(gin_fr, "FRENERTRJ");
	if(end_of_file){
	  if(it_fr!=to_fr){
	    gin_fr.close();
	    gin_fr.open(it_fr->second.c_str());
	    ++it_fr;
	    //try again...
	    etrj.read_frame(gin_fr, "FRENERTRJ");
	  }
	  else
	    break;
	}
      }

      int eblockindex = etrj.return_blockindex("PRECALCLAM");
      int feblockindex = etrj.return_blockindex("FREEPRECALCLAM");
      vector<vector<double> > se;
      vector<vector<double> > sfe;
      vector<double> seRef;
      vector<double> sfeRef;

      // read data
      if(!single_plam){
        se = etrj.return_block(eblockindex);
        sfe = etrj.return_block(feblockindex);
      }
      else {
        se.resize(1);
        sfe.resize(1);
        se[0] = etrj.return_line(eblockindex,pindex);
        sfe[0] = etrj.return_line(feblockindex,pindex);
      }

      // read ref data
      seRef = etrj.return_line(eblockindex,sindex);
      sfeRef = etrj.return_line(feblockindex,sindex);

      // se[l][0] : A_e_lj
      // se[l][1] : B_e_lj
      // se[l][2] : A_e_crf
      // se[l][3] : B_e_crf
      // sfe[l][0] : A_de_lj
      // sfe[l][1] : B_de_lj
      // sfe[l][2] : A_de_crf
      // sfe[l][3] : B_de_crf

      double dlamlj_dlam = 1.0;
      double dlamslj_dlam = 1.0;
      double dlamcrf_dlam = 1.0;
      double dlamscrf_dlam = 1.0;

      // calculate the reference energies
      double pow1_sl_ns = pow((1-slam),NLAMs);
      double pow_sl_ns = pow(slam,NLAMs);
      double Elj_s = pow1_sl_ns * seRef[0] + pow_sl_ns * seRef[1];
      double Ecrf_s = pow1_sl_ns * seRef[2] + pow_sl_ns * seRef[3];
      double p_lam;

      // for all plam, calculate the energies
      for(int p=0; p<nr_plam; p++){
        if(single_plam) 
          p_lam = plam; 
        else
          p_lam = p *lam_step + minlam;

        // precalculate some constants
        double pow1_pl_np = pow((1-p_lam),NLAMp);
        double pow_pl_np = pow(p_lam,NLAMp);
        double np_pow1_pl_np_1 = NLAMp * pow((1-p_lam),(NLAMp-1));
        double np_pow_pl_np_1 = NLAMp * pow(p_lam,(NLAMp-1));

        // calculate reference energies
        double Elj = pow1_pl_np * se[p][0] + pow_pl_np * se[p][1];
        double Ecrf = pow1_pl_np * se[p][2] + pow_pl_np * se[p][3];

        // and their derivatives
        double dElj = (-np_pow1_pl_np_1 * se[p][0] + np_pow_pl_np_1 * se[p][1] ) * dlamlj_dlam + (pow1_pl_np * sfe[p][0] - pow_pl_np * sfe[p][1] ) * dlamslj_dlam;
        double dEcrf = (-np_pow1_pl_np_1 * se[p][2] + np_pow_pl_np_1 * se[p][3] ) * dlamcrf_dlam + (pow1_pl_np * sfe[p][2] - pow_pl_np * sfe[p][3] ) * dlamscrf_dlam;

        // prepare for reweighting
        double X = dElj+dEcrf;
        double diff = -(Elj - Elj_s + Ecrf - Ecrf_s)/(gmath::physConst.get_boltzmann() * temp);
        x[p].addval(X);
        double expdiff = exp(diff);
        vyvr[p].addval(diff);
        expvyvr[p].addval(expdiff);
        Xexpvyvr[p].addval(X * expdiff);
      }
    }


    // for each plam
    for(int p=0; p<nr_plam; p++){
      double p_lam;
      if(plam)
        p_lam = plam;
      else
        p_lam = p *lam_step + minlam;

      if(printX && (p_lam == print_plam)){
        ostringstream os;
        os << "X.dat";
        ofstream fout(os.str().c_str());
        fout.precision(9); //set precision of numbers going to ofstream
        fout << "#"
             << setw(14) << "X"
             << endl;
        for(int i=0; i< x[p].n(); i++){
          fout << setw(15) << x[p].val(i)
               << endl;
        }
        fout.close();

        // also print expvyvr
        ostringstream os1;
        os1 << "expvyvr.dat";
        ofstream fout1(os1.str().c_str());
        fout1.precision(9);
        fout1 << "#"
             << setw(14) << "expvyvr"
             << endl;
        for(int i=0; i< expvyvr[p].n(); i++){
          fout1 << setw(15) << expvyvr[p].val(i)
               << endl;
        }
        fout1.close();
      }

      // calculate the final value
      int sign = 0;
      double lnXexpave = gmath::Stat<double>::lnXexpave(x[p], vyvr[p], sign);
      double lnexpave = vyvr[p].lnexpave();
      double final = exp(lnXexpave - lnexpave) * sign;

      // calculate statistical uncertainty if required
      if(no_error){
        cout << p_lam << setw(15) << final << endl;
      }
      else{
        double dii = exp(2 * lnXexpave);
        double djj = exp(2 * lnexpave);
        double dji = exp(lnXexpave + lnexpave);
       
        double n = 1.0 / x[p].n();
        double var_ii = gmath::Stat<double>::covariance(Xexpvyvr[p], Xexpvyvr[p]);
        double si_ii = gmath::Stat<double>::stat_ineff(Xexpvyvr[p], Xexpvyvr[p]);
        double d2i = var_ii * si_ii * n;

        double error;
        // if plam equals slam, we only calculate the error estimate
        // for Xexpvyvr
        if( slam == p_lam){
          error = sqrt( d2i );
        }
        else{
          double var_jj = gmath::Stat<double>::covariance(expvyvr[p], expvyvr[p]);
          double si_jj = gmath::Stat<double>::stat_ineff(expvyvr[p], expvyvr[p]);
          double d2j = var_jj * si_jj * n;
          double var_ji = gmath::Stat<double>::covariance(Xexpvyvr[p], expvyvr[p]);
          double si_ji = gmath::Stat<double>::stat_ineff(Xexpvyvr[p], expvyvr[p]);
          double d2ji = var_ji * si_ji * n;
          error = abs(final) * sqrt( d2i/dii + d2j/djj - 2*d2ji/dji );
//          if( isnan(error) ){
//            error = abs(final) * sqrt(abs(d2i/dii + d2j/djj - 2*d2ji/dji));
//            if( isnan(error) ){
//              error = 100.;
//            }
//          }
        }
        cout << p_lam << setw(15) << final << setw(15) 
             << error << endl;

//        double dii = exp(2 * lnXexpave);
//        double djj = exp(2 * lnexpave);
//        double dji = exp(lnXexpave + lnexpave);
       
//        double n = 1.0 / x[p].n();
//        double var_ii = gmath::Stat<double>::covariance(Xexpvyvr[p], Xexpvyvr[p]);
//        double si_ii = gmath::Stat<double>::stat_ineff(Xexpvyvr[p], Xexpvyvr[p]);
//        double d2i = var_ii * si_ii * n;
//        double var_jj = gmath::Stat<double>::covariance(expvyvr[p], expvyvr[p]);
//        double si_jj = gmath::Stat<double>::stat_ineff(expvyvr[p], expvyvr[p]);
//        double d2j = var_jj * si_jj * n;
//        double var_ji = gmath::Stat<double>::covariance(Xexpvyvr[p], expvyvr[p]);
//        double si_ji = gmath::Stat<double>::stat_ineff(Xexpvyvr[p], expvyvr[p]);
//        double d2ji = var_ji * si_ji * n;
//        double error = sqrt( d2i/dii + d2j/djj - 2*d2ji/dji );
       
//        cout << p_lam << setw(15) << final << setw(15) 
//             << error << endl;
      }
    }
    
  }

  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


void set_standards(utils::EnergyTraj &e)
{  
  e.addConstant("BOLTZ", gmath::physConst.get_boltzmann());
}

void read_library(string name, utils::EnergyTraj& e)
{
  Ginstream gin;
  
  try{
    
    gin.open(name);
  }
  catch (gromos::Exception ex){
      throw gromos::Exception("read_library", "failed to open library file "
			      +name);
  }
  while(true){
    
    vector<string> buffer;
    gin.getblock(buffer);
    if(gin.stream().eof()) break;
    if(buffer[buffer.size()-1].find("END")!=0)
      throw gromos::Exception("ene_ana", "Library file " + gin.name() +
			      " is corrupted. No END in "+buffer[0]+
			      " block. Got\n"
			      + buffer[buffer.size()-1]);
    string sdum;
    
    if(buffer[0]=="ENERTRJ" || buffer[0]=="FRENERTRJ"){
      for(unsigned int i=1; i< buffer.size()-1; i++){
	e.addBlock(buffer[i], buffer[0]);
	
      }
    }
    

    vector<string> data;
    if(buffer[0]=="VARIABLES"){
      
      set_standards(e);
      
      string bufferstring;
      
      gio::concatenate(buffer.begin()+1, buffer.end(), bufferstring);
      
      istringstream iss(bufferstring);

      // i am aware of the fact that END will also be stored in data.
      // This is used in parsing later on
      while(sdum!="END"){
	iss >> sdum;
	data.push_back(sdum);
      }
      
      // now search for the first appearance of "="
      for(unsigned int i=0; i< data.size(); i++){
	if(data[i]=="="){
	  
	  // search for either the next appearance or the end
	  unsigned int to=i+1;
	  for(; to < data.size(); to++) if(data[to]=="=") break;
	  
	  // parse the expression part into an ostringstream
	  ostringstream os;
	  for(unsigned int j=i+1; j< to-1; j++) os << " " << data[j]; 
	  e.addKnown(data[i-1], os.str());
	}
      }
    }
  }
}
