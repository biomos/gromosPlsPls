// dg_ener.cc This program calculates a free energy difference by reading
//            in the output files of ener, for state A and B

#include "../src/args/Arguments.h"
#include "../src/gmath/physics.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <sstream>

using namespace args;
using namespace std;

int main(int argc, char **argv){

  char *knowns[] = {"temp", "stateA", "stateB"};
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@temp <temperature for perturbation>\n";
  usage += "\t@stateA <energy file for state A>\n";
  usage += "\t@stateB <energy file for state B>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // set some values
    args.check("temp", 1);
    double temp=0;
    
    {
      Arguments::const_iterator iter=args.lower_bound("temp");
      if(iter!=args.upper_bound("temp"))
        temp=atof(iter->second.c_str());
    }
    ifstream stateA, stateB;

    //open files
    {
      Arguments::const_iterator iter=args.lower_bound("stateA");
      if(iter!=args.upper_bound("stateA")){
        stateA.open((iter->second).c_str());
        if(!stateA)
          throw gromos::Exception("dg_ener", "could not open energy file 
for state A\n"); 
      }
      else
        throw gromos::Exception("dg_ener", "energy file for state A missing\n");
      iter=args.lower_bound("stateB");
      if(iter!=args.upper_bound("stateB")){
	stateB.open((iter->second).c_str());
        if(!stateB)
          throw gromos::Exception("dg_ener", "could not open energy file 
for state B\n"); 
      }
      else
	throw gromos::Exception("dg_ener", "energy file for state B missing\n");
    }
    cout.precision(12);
    
    //print title
    cout << "# Time"
	 << setw(12) << "DE_tot"
	 << setw(12) << "probability"
	 << setw(12) << "DG_AB" 
	 << endl;
    
    string sdum;
    double timeA, timeB;
    
    bool errorA = false, errorB = false;
    bool eofA = false, eofB = false;
    bool timeWarning = false;

    vector<double> delta_v;
    vector<double> t;
    double ave_v=100000.0;

    double covA, nbA, totA, covB, nbB, totB;

    while(true){

      // read from state A
      while(true){
	std::getline(stateA, sdum);
	if (stateA.eof()){
	  eofA = true;
	  break;
	}
	std::string::size_type it = sdum.find('#');
	if (it != std::string::npos)
	  sdum = sdum.substr(0, it);
	if (sdum.find_first_not_of(" \t") == std::string::npos)
	  continue;
	std::istringstream is(sdum);
	if (!(is >> timeA >> covA >> nbA >> totA))
	  errorA = true;
	break;
      }
      // read from state B
      while(true){
	std::getline(stateB, sdum);
	if (stateB.eof()){
	  eofB = true;
	  break;
	}
	std::string::size_type it = sdum.find('#');
	if (it != std::string::npos)
	  sdum = sdum.substr(0, it);
	if (sdum.find_first_not_of(" \t") == std::string::npos)
	  continue;
	std::istringstream is(sdum);
	if (!(is >> timeB >> covB >> nbB >> totB))
	  errorB = true;
	break;
      }
      // check eof / error
      if (eofA && eofB) break;
      if (eofA || eofB || errorA || errorB)
	throw gromos::Exception("dg_ener", "Error while reading file for state A or state B: check if number of lines or columns are identical");
      if (timeA != timeB)
	timeWarning = true;
      
      // now we have two states read
      const double dh = totB - totA;
      delta_v.push_back(dh);
      
      if (dh < ave_v) ave_v = dh;
      t.push_back(timeA);
    }

    double sum=0, p, dg=0.0, ave;
    
    // now loop over all values that were read in 
    for(unsigned int i=0; i<delta_v.size(); i++){

      const double dh=delta_v[i] - ave_v;
      p = exp(-dh/(BOLTZ * temp));
      sum += p;
      ave = sum / i;
      dg = ave_v - BOLTZ * temp * log(ave);
      
      cout.precision(5);
      cout.setf(ios::right, ios::adjustfield);
      cout << setw(6) << t[i]
           << setw(12) << delta_v[i]
           << setw(12) << p
           << setw(12) << dg
           << endl;
    }

    cout << "# substracted value " << ave_v << endl;
    cout << "# final result: " << dg << endl;

    stateA.close();
    stateB.close();

    if (timeWarning){
      cout << "#\n# WARNING: time was not equal in state A and state B!" << endl;
      cerr << "\nWARNING: time was not equal in state A and state B!" << endl;
    }
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    return 1;
  }
  return 0;
}
