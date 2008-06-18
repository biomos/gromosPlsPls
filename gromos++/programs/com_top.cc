/**
 * @file com_top.cc
 * Combine molecular topology files into one
 */

/**
 * @page programs Program Documentation
 *
 * @anchor com_top
 * @section com_top Combine molecular topology files into one
 * @author @ref co
 * @date 7-6-07
 *
 * To generate molecular topology files for the use in simulations of e.g. 
 * (macro)molecular complexes, or mixtures containing several solutes and/or 
 * (co)solvents, it is usually convenient to merge existing molecular topology
 * files. Program com top combines multiple topologies into one new topology.
 * 
 * The user has to specify which molecular topologies to be merged (use prefix
 * 'n:' before the file name to repeat one topology n times), and from  which 
 * file the force field parameters and the solvent have to be taken. The
 * resulting molecular topology file is written out to the standard output. 
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology files&gt; </td></tr>
 * <tr><td> \@param</td><td>&lt;index number of molecular topology file to take parameters from&gt; </td></tr>
 * <tr><td> \@solv</td><td>&lt;index number of molecular topology file to take solvent from&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  com_top
    @topo    ex.top 7:cl.top
    @param   1
    @solv    1
 @endverbatim
 *
 * <hr>
 */


#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "param", "solv"};
  int nknowns = 3;

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology files>\n";
  usage += "\t@param <index number of molecular topology file to take parameters from>\n";
  usage += "\t@solv  <index number of molecular topology file to take solvent from>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    // set some values
    int parnum=1;
    if(args.count("param")>0) parnum=atoi(args["param"].c_str());
    int solnum=1;
    if(args.count("solv")>0)  solnum=atoi(args["solv"].c_str());
    if(args.count("topo")<=0)
      throw gromos::Exception("com_top", "needs at least one topology\n"+usage);
    
    System sys;

    //ugly solution to the not yet implemented '=' for force fields
    string paramname, toponame, s;
    std::string::size_type s_it;
    int repeat = 1, topnum=0, oldtopnum=-1;
    
    ostringstream title;
    title << "COM_TOP: Combined topology using:\n";
    
   
    int totNumAt=0;    
    for(Arguments::const_iterator iter=args.lower_bound("topo"),
         to=args.upper_bound("topo"); iter!=to; ++iter){
      oldtopnum=topnum+1;
      
      s=iter->second;
      s_it=s.find(':');
      
      if(s_it == string::npos){
	toponame=s;
	repeat=1;
      }
      else{
	toponame=s.substr(s_it+1,s.size());
	repeat=atoi(s.substr(0,s_it).c_str());
      }
      topnum+=repeat;
      
      // read topology
      InTopology it(toponame);
      for(int i=0; i<repeat; i++){
   
    // Directly add pressure and temperature groups
    for(int j=0; j<it.system().numTemperatureGroups(); j++){
      sys.addTemperatureGroup(it.system().temperatureGroup(j)+totNumAt);
    }
    for(int j=0; j<it.system().numPressureGroups(); j++){
      sys.addPressureGroup(it.system().pressureGroup(j)+totNumAt);
    }

    // Add molecules and count new number of atoms in sys
	for(int j=0;j<it.system().numMolecules();j++){
	  sys.addMolecule(it.system().mol(j));
      totNumAt+=it.system().mol(j).numAtoms();
    }
      }

      if(solnum <= topnum && solnum >= oldtopnum)
	sys.addSolvent(it.system().sol(0));
      if(parnum <= topnum && parnum >= oldtopnum)
        paramname=toponame;
      if(topnum!=oldtopnum)
	title << setw(4) << oldtopnum << " .. " << setw(4) << topnum;
      else
	title << setw(12) << topnum;
      
      title << " : " << toponame << endl;
    }
    InTopology it(paramname);
    title << "Parameters from " << parnum 
          << ", solvent from " << solnum;
    
    OutTopology ot(cout);
    
    ot.setTitle(title.str());
    ot.write(sys,it.forceField());
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




