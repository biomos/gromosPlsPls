/**
 * @file ene_ana.cc
 * extracts time series from (energy) trajectory files
 */

/**
 * @page programs Program Documentation
 *
 * @anchor ene_ana
 * @section ene_ana analyse (energy) trajectories
 * @author @ref mc @ref co
 * @date 26. 7. 2006
 *
 * GROMOS can write energies, free-energy derivatives and block averages of 
 * these to separate trajectory files for later analysis. Program ene_ana 
 * extracts individual values from such files and can perform simple 
 * mathematical operations on them. The format for (free) energy trajectory 
 * files as written by promd and md are known to the program. In addition, the
 * user can define custom made formats of any trajectory file that comes in a
 * block-format through a library file. ene_ana is able to read and interpret 
 * series of two types of such files simultaneously, typically referred to as 
 * the "energy file" and the "free energy file".
 * 
 * Using the same library file one can define properties to be calculated from 
 * the values that are listed in them. For the selected properties, ene_ana 
 * will calculate the time series, averages, root-mean-square fluctuations and 
 * a statistical error estimate. The error estimate is calculated from block 
 * averages of different sizes, as described in Allen and Tildesley: "Computer 
 * Simulation of Liquids", 1987. The time for the time series is taken from the
 * trajectory files, unless a different time interval between blocks is 
 * specified through an input parameter. If a topology is supplied, the ene_ana 
 * uses this to define the total solute mass (MASS) and the total number of 
 * solute molecules (NUMMOL).
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@en_files</td><td>&lt;energy files&gt; (and/or) </td></tr>
 * <tr><td> \@fr_files</td><td>&lt;free energy files&gt; </td></tr>
 * <tr><td> \@prop</td><td>&lt;@ref PropertySpecifier "properties" to monitor&gt; </td></tr>
 * <tr><td> \@library</td><td>&lt;library for property names&gt; [print] </td></tr>
 * <tr><td> [\@topo</td><td>&lt;molecular topology file&gt; (for MASS and NUMMOL)] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt; (overwrites TIME in the trajectory files)] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ene_ana
    @topo       ex.top
    @en_files   ex.tre
    @prop       densit
    @library    ene_ana.lib

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

void print(gmath::Stat<double> &p, string s, vector<double> & time);
void set_standards(utils::EnergyTraj &e);
void read_library(string name, utils::EnergyTraj& e);

int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "time" << "en_files" << "fr_files" << "prop" << "library";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@en_files    <energy files> (and/or)\n";
  usage += "\t@fr_files    <free energy files>\n";
  usage += "\t@prop        <properties to monitor>\n";
  usage += "\t@library     <library for property names> [print]\n";
  usage += "\t[@topo       <molecular topology file> (for MASS and NUMMOL)]\n";
  usage += "\t[@time       <t and dt> (overwrites TIME in the trajectory files)]\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    // get simulation time either from the user or from the files
    bool usertime=false;
    vector<double> time;
    double t0=0, dt=1; 
    {
      Arguments::const_iterator iter=args.lower_bound("time");
      if(iter!=args.upper_bound("time")){
        t0=atof(iter->second.c_str());
        ++iter;
      }
      if(iter!=args.upper_bound("time")){
        dt=atof(iter->second.c_str());
	usertime=true;
        // a bit ugly: as the time is increased by dt before the first printout: reduce t0
	t0 -= dt;
      }
    }

    // check whether we are doing anything
    if(args.count("en_files")<=0 && args.count("fr_files")<=0)
      throw gromos::Exception("ene_ana", "no data specified:\n"+usage);
    if(args.count("prop") <=0)
      throw gromos::Exception("ene_ana", "no properties to follow:\n"+usage);
    
    // NEW: require a library 
    if(args.count("library") <=0)
      throw gromos::Exception("ene_ana", "no library file specified:\n"+usage);
    
    // read a library file?
    string library="";
    int print_library=0;
    {
      Arguments::const_iterator iter=args.lower_bound("library"), 
	to=args.upper_bound("library");
      if(iter!=to){
	library=iter->second;
	++iter;
      }
      if(iter!=to && iter->second == "print") print_library=1;
    }
    
    // which properties do we follow
    vector<string> prop;
    int num_prop=0;
    {
      Arguments::const_iterator iter=args.lower_bound("prop"), 
	to=args.upper_bound("prop");
      while(iter!=to){
	prop.push_back(iter->second);
	++iter;
      }
      num_prop=prop.size();
    }

    // define an energy trajectory
    utils::EnergyTraj etrj;

    // read topology for the mass and the number of molecules
    double mass=0;
    if(args.count("topo")>0){
      
      InTopology it(args["topo"]);
      System sys(it.system());
      etrj.addConstant("NUMMOL", sys.numMolecules());
      // calculate total mass of the system
      for(int m=0; m<sys.numMolecules(); m++){
	for(int a=0; a< sys.mol(m).topology().numAtoms(); a++){
	  mass += sys.mol(m).topology().atom(a).mass();
	}
      }
    }
    etrj.addConstant("MASS",mass);
    
    // learn about the variable names how they map to the elements
    read_library(library, etrj);
    
    if(print_library) etrj.write_map();

    // prepare for the statistical information
    vector<gmath::Stat<double> > s(num_prop);

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
    
    bool version_checked = false;
    while (true) {
      
      // version number
      if (!version_checked) {
        version_checked = true;
        if (etrj.has_version()) {
          if (do_energy_files) {
            if (gin_en.has_version()) {
              if (!etrj.version_match(gin_en.version())) {
                cerr << "WARNING: Version number mismatch!\n"
                     << "         Library " << library << " version: "
                     << etrj.get_version() << std::endl
                     << "         Energy Trajectory " << gin_en.name() << " version: "
                     << gin_en.version() << std::endl;
              } else {
                cerr << "MESSAGE: Version number check successful!\n"
                     << "         Library " << library << " version: "
                     << etrj.get_version() << std::endl
                     << "         Energy Trajectory " << gin_en.name() << " version: "
                     << gin_en.version() << std::endl;
              }
            } else {
              cerr << "WARNING: Version number missing!\n"
                   << "         Library " << library << " version: "
                   << etrj.get_version() << std::endl
                   << "         Energy Trajectory " << gin_en.name() << " version: "
                   << "NONE" << std::endl;
            }
          }
          if (do_free_energy_files) {
            if (gin_fr.has_version()) {
              if (!etrj.version_match(gin_fr.version())) {
                cerr << "WARNING: Version number mismatch!\n"
                     << "         Library " << library << " version: "
                     << etrj.get_version() << std::endl
                     << "         Energy Trajectory " << gin_fr.name() << " version: "
                     << gin_fr.version() << std::endl;
              } else {
                cerr << "MESSAGE: Version number check successful!\n"
                     << "         Library " << library << " version: "
                     << etrj.get_version() << std::endl
                     << "         Energy Trajectory " << gin_fr.name() << " version: "
                     << gin_fr.version() << std::endl;
              }
            } else {
              cerr << "WARNING: Version number missing!\n"
                   << "         Library " << library << " version: "
                   << etrj.get_version() << std::endl
                   << "         Energy Trajectory " << gin_fr.name() << " version: "
                   << "NONE" << std::endl;
            }
          }
        } else {
          if (do_energy_files) {
            if (gin_en.has_version()) {
              cerr << "WARNING: Version number missing!\n"
                   << "         Library " << library << " version: "
                   << "NONE" << std::endl
                   << "         Energy Trajectory " << gin_en.name() << " version: "
                   << gin_en.version() << std::endl;
            } else {
              cerr << "WARNING: Version number missing!\n"
                   << "         Library " << library << " version: "
                   << "NONE" << std::endl
                   << "         Energy Trajectory " << gin_en.name() << " version: "
                   << "NONE" << std::endl;
            }
          }
          if (do_free_energy_files) {
            if (gin_fr.has_version()) {
              cerr << "WARNING: Version number missing!\n"
                   << "         Library " << library << " version: "
                   << "NONE" << std::endl
                   << "         Energy Trajectory " << gin_fr.name() << " version: "
                   << gin_fr.version() << std::endl;
            } else {
              cerr << "WARNING: Version number missing!\n"
                   << "         Library " << library << " version: "
                   << "NONE" << std::endl
                   << "         Energy Trajectory " << gin_fr.name() << " version: "
                   << "NONE" << std::endl;
            }
          }
        }
      }      
      
      // read the numbers into the energy trajectory
      if(do_energy_files){
	int end_of_file=etrj.read_frame(gin_en, "ENERTRJ");
	if(end_of_file){
	  if(it_en!=to_en){
	    gin_en.close();
	    gin_en.open(it_en->second.c_str());
	    ++it_en;
	    //try again...
	    etrj.read_frame(gin_en, "ENERTRJ");
        version_checked = false;
	  } else {
	    if(do_free_energy_files){
	      // check if also the free energy traj is finished
	      int end_of_file=etrj.read_frame(gin_fr, "FRENERTRJ");
	      if(!end_of_file){
	         std::cerr << "# WARNING: frames left over in free energy trajectory,\n"
	                 << "#   they will be disregarded (ene_ana processes data\n"
	                 << "#   frame-by-frame in parallel for energy and free energy trajectories).\n"
	                 << "#   Probably you used different timesteps for trg and tre writeout.\n";
	      }
	    }
	    break;
	  }
	}
      }
      if(do_free_energy_files){
	int end_of_file=etrj.read_frame(gin_fr, "FRENERTRJ");
	if(end_of_file){
	  if(it_fr!=to_fr){
	    gin_fr.close();
	    gin_fr.open(it_fr->second.c_str());
	    ++it_fr;
	    //try again...
	    etrj.read_frame(gin_fr, "FRENERTRJ");
        version_checked = false;
	  } else {
        if(do_energy_files){
          // we are only ending up here if the energy trajectory has less frames
	      std::cerr << "# WARNING: frames left over in energy trajectory,\n"
	                 << "#   they will be disregarded (ene_ana processes data\n"
	                 << "#   frame-by-frame in parallel for energy and free energy trajectories).\n"
	                 << "#   Probably you used different timesteps for trg and tre writeout.\n";
	    }
	    break;
	  }
	}
      }
      
      
      // calculate and store the necessary number in the stat-classes
      for(int i=0; i<num_prop; i++)
	s[i].addval(etrj[prop[i]]);
      if(usertime)
	t0+=dt;
      else      
	t0=etrj["TIME[2]"];
      time.push_back(t0);
    }
    //print out the statistical information
    cout << setw(10) << "property"
     << " "
	 << setw(15) << "average"
	 << " "
	 << setw(15) << "rmsd"
	 << " "
	 << setw(15) << "error est."
	 << endl;
    for(int i=0; i<num_prop; i++)
      print(s[i], prop[i], time);
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



void print(gmath::Stat<double> &p, string s, vector<double>& time)
{
  if(p.n()!=int(time.size())) 
    throw gromos::Exception("ene_ana","number of time steps read is not equal"
			    " to number of data points to print out");
  
  ostringstream os;
  os << s << ".dat";
  ofstream fout(os.str().c_str());
  fout.precision(9); //set precision of numbers going to ofstream
  fout << "#"
       << setw(14) << "time"
       << " "
       << setw(15) << s
       << endl;
  for(int i=0; i< p.n(); i++){
    fout << setw(15) << time[i]
     << " "
	 << setw(15) << p.val(i)
	 << endl;
  }
  fout.close();
// and print the averages etc to cout
  cout.precision(9); // set precision of number going to cout
  cout << setw(10) << s
       << " " //put an extra space, that will always pe printed, even if setw space does not suffice, to prevent merging of columns
       << setw(15) << p.ave()
       << " "
       << setw(15) << p.rmsd()
       << " "
       << setw(15) << p.ee()
       << endl;
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
    
    if (buffer[0] == "ENEVERSION") {
      // This is a ENEVERSION block that was not put directly after the
      // TITLE block. This is no big deal for the library, it would
      // however lead to an error in a trajectory.
      // Shall we disallow this?
      
      // Just to make sure there are not two blocks
      if (gin.has_version() || e.has_version()) {
        throw gromos::Exception("ene_ana", "Library file " + gin.name() +
            " is corrupted. Two ENEVERSION blocks found.");
      }
      string ver;
      gio::concatenate(buffer.begin() + 1, buffer.end()-1, ver);
      ver.erase( std::remove_if( ver.begin(), ver.end(), ::isspace ), ver.end() );
      e.set_version(ver);
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
  if (gin.has_version()) {
    e.set_version(gin.version());
  }
}
