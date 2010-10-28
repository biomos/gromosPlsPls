/**
 * @file close_pair.cc
 * find atom pairs of two molecules that are close to each other
 */

/**
 * @page programs Program Documentation
 *
 * @anchor close_pair
 * @section close_pair
 * @author @ref dw
 * @date 01. 07. 2010
 *
 * Program close_pair finds the closest atom pairs of two molecules.
 * The searching for molecule i loops over all solute atoms of the other molecules,
 * and if the distance between the i:j and k:l is shorter than a value (dc default 0.3 nm),
 * the searching for the atom pairs between i - k is stopped. Periodicity is
 * also considered if the pbc is defined.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@atomsrmsf</td><td>&lt;@ref AtomSpecifier "atoms" to consider for rmsf&gt; </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;@ref AtomSpecifier "atoms" to consider for fit&gt;] </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinates(if absent, the first frame of \@traj is reference)&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
  rsmf
    @topo       ex.top
    @pbc        r
    @groupA   1
    @groupB   2
    @traj       ex.tr
@endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <iostream>
#include <time.h>

#include "../src/args/Arguments.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/LJExceptionType.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/groTime.h"
#include "../src/gcore/Box.h"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;



int main(int argc, char **argv){

  Argument_List knowns; 
  knowns <<"topo" << "traj" << "groupA" << "groupB" << "pbc" << "time" << "dist" << "debug";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc         <boundary type>\n";
  usage += "\t@groupA      <molecule group to be gathered>\n";
  usage += "\t@groupB      <molecule group to be reference for gathering of groupA>]\n";
  usage += "\t[@dist       <lower limit of distance dc. default 0.2 nm>]\n";
  usage += "\t[@time       <t0 and dt>]\n";
  usage += "\t@traj        <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);


  try{
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    //System refSys(it.system());
    
    // read reference coordinates...
    InG96 ic;

            //get time
    utils::Time time(args);
    //Life solvlife(args);

    /*
    try{
      args.check("ref",1);
      ic.open(args["ref"]);
    }

    catch(const Arguments::Exception &){
      args.check("traj",1);
      ic.open(args["traj"]);
    }
    ic >> refSys;
    ic.close();
     */
    double dist = 0.30;
    double dc;
    try{
      args.check("dist", 1);
      dist = atof(args["dist"].c_str());
    }
    catch (const gromos::Exception &e){}
    // we will calc only the square of distance
    dc = dist * dist;

    // get simulation time
    /*
    double time0=0, dt=1;
    {
      Arguments::const_iterator iter=args.lower_bound("time");
      if(iter!=args.upper_bound("time")){
        time0=atof(iter->second.c_str());
        ++iter;
      }
      if(iter!=args.upper_bound("time"))
        dt=atof(iter->second.c_str());
    }
     */

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    
    // System for calculation
    //System sys(refSys);
    //get the atoms for the rmsf
    AtomSpecifier groupA(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("groupA");
      Arguments::const_iterator to = args.upper_bound("groupA");
      
      for(;iter!=to;iter++){
	string spec=iter->second.c_str();
        //if(solvatoms.name(spec) == "OW")
        groupA.addSpecifier(spec);
      }
    }

    AtomSpecifier groupB(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("groupB");
      Arguments::const_iterator to = args.upper_bound("groupB");

      for(;iter!=to;iter++){
	string spec=iter->second.c_str();
        //if(solvatoms.name(spec) == "OW")
        groupB.addSpecifier(spec);
      }
    }
    int debug=0;
    if(args.count("debug")>0){
        debug=atoi(args["debug"].c_str());
    }
    if(debug==1){
        for(int i=0;i<groupA.size();++i){
            cout << "# groupA mol " << groupA.mol(i) << " size " << groupA.size() << endl;
        }
        for(int j=0;j<groupB.size();++j){
            cout << "# groupB mol " << groupB.mol(j) << " size " << groupB.size() << endl;
        }
    }

    if (groupA.size() == 0)
        throw gromos::Exception("groupA",
                "No atoms specified for groupA");

    groupA.sort();
    groupB.sort();

    double boxsquare=sys.box().X();
    boxsquare *=boxsquare/4.0;

    cout << setw(10) << "#    Frame" << setw(8) << "molA"
            << setw(8) << "resA" << setw(6) << "atomA"
            << setw(8) << "molB" << setw(8) << "resB" << setw(6) << "atomB"
            << setw(10) << "atomA" << setw(8) << "atomB"
            << setw(12) << "dist [nm]" << endl;

    int numFrames = 0;
    
    //vector<Vec> apos;
    //Vec apos(0.,0.,0.);
    vector<double> apos2;
    vector<int> aresid, aatmid;
    vector<string> aatom;
    vector<int> bmolid, bresid, batmid;
    vector<string> batom;

    //dimension= sys.numMolecules() * sys.numMolecules();
    apos2.resize(sys.numMolecules(), 0.0); // shortest dist
    aresid.resize(sys.numMolecules(),0);  // molA:ResID
    aatmid.resize(sys.numMolecules(),0);  // molA:AtomID
    aatom.resize(sys.numMolecules(),"");  // molA:AtomName

    bmolid.resize(sys.numMolecules(),0);
    bresid.resize(sys.numMolecules(),0);
    batmid.resize(sys.numMolecules(),0);
    batom.resize(sys.numMolecules(),"");

    // vector to save temporary values
    vector< vector<double> > buff_apos2(sys.numMolecules(),apos2);
    vector< vector<int> > buff_aresid(sys.numMolecules(),aresid), buff_aatmid(sys.numMolecules(),aatmid);
    vector< vector<string> > buff_aatom(sys.numMolecules(),aatom);
    vector< vector<int> > buff_bmolid(sys.numMolecules(),bmolid), buff_bresid(sys.numMolecules(),bresid), buff_batmid(sys.numMolecules(),batmid);
    vector< vector<string> > buff_batom(sys.numMolecules(),batom);

    double temp; // save the dist^2 and then compare with apos2, save the shorter one

    double clo;

    //double diagonal=sys.box().X() * sys.box().X() + sys.box().Y() * sys.box().Y() + sys.box().Z() * sys.box().Z();

    //loop over trajectory
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      // loop over all frames
      while(!ic.eof()){
        // read frame
        ic.select("ALL");
	//ic >> sys >> time;
        ic >> sys ;

        if(debug){
            clo = clock();
        }

        // initalize after it is read
        //if (numFrames == 0) {
          
        int pass=1000;

	for(int i=0; i< groupA.size(); ++i){
            //solvindex = solvatoms.gromosAtom(i);
          // skip hydrogen atoms
            int m = groupA.mol(i);

            if(debug>=3)
                cout << "# search for atom  " << i << "   " << m << ":" << groupA.resnum(i) << " " << groupA.name(i) << endl;

            if(groupA.resname(i)=="")
                      throw gromos::Exception("groupA",
                              "atom name not specified properly.");
            //if(debug)
            //    cout << "# working on groupA mol " << m << "\t atom " << i << " name " << groupA.name(i) << endl;
            if(m!=pass && m>0)
            if(groupA.mass(i) != 1.008){
                for(int j=0; j<groupB.size(); ++j){
                    if(debug>=3){
                        cout << "# groupB size " << groupB.size() << endl;
                        cout << "# molA " << groupA.mol(i)
                                << " atom " << i
                                << " name " << groupA.name(i)
                                << " mass " << groupA.mass(i)
                                << "\t molB " << groupB.mol(j)
                                << " atom " << j
                                << " name " << groupB.name(j)
                                << " mass " << groupB.mass(j) << endl;
                    cout << "# groupB.mol(j) " << groupB.mol(j) << " groupA.mol(i) " << groupA.mol(i) << endl;
                    }
                    if(groupB.name(j) == "")
                        throw gromos::Exception("groupB",
                              "atom name not specified properly.");
                    //cout << "# groupB.mol(j) " << groupB.mol(j) << " groupA.mol(i) " << groupA.mol(i) << endl;
                    if(groupB.mol(j)<groupA.mol(i))
                        if(groupB.mass(j)!=1.008){
                        //apos2[i] = (groupA.pos(i) - groupB.pos(j)).abs2();
                            //temp=(groupA.pos(i) - groupB.pos(j)).abs2();
                            Vec ref=pbc->nearestImage(groupA.pos(i),groupB.pos(j),sys.box());
                            temp=(groupA.pos(i) - ref).abs2();
                            if(debug)
                                cout << "# done for nearest image " << j << endl;
                            if(debug){
                                double distance=(groupA.pos(i) - groupB.pos(j)).abs2();
                                cout << m << " : " << i << "\t" << groupB.mol(j) << " : " << j
                                    << "\t" << temp << "\t" << apos2[m] << "\t" << distance << endl;
                            }

                            if(debug>=3)
                                cout << "# done for pair: " << m << ":" << i << " " << groupA.name(i) << "\t"
                                    << groupB.mol(j) << ":" << j << " " << groupB.name(j) << "\t" << temp << endl;

                            if(buff_apos2[m][groupB.mol(j)]==0 && temp!=0){
                                buff_apos2[m][groupB.mol(j)] = temp;
                                buff_aresid[m][groupB.mol(j)] =groupA.resnum(i);
                                buff_aatmid[m][groupB.mol(j)] =groupA.atom(i);
                                buff_aatom[m][groupB.mol(j)] =groupA.name(i);

                                buff_bmolid[m][groupB.mol(j)] =groupB.mol(j);
                                buff_bresid[m][groupB.mol(j)] =groupB.resnum(j);
                                buff_batmid[m][groupB.mol(j)] =groupB.atom(j);
                                buff_batom[m][groupB.mol(j)] =groupB.name(j);
                                if(debug)
                                    cout << "buff term m " << m << " molB " << groupB.mol(j) << " initialized " << endl;
                            }else{
                              if(buff_apos2[m][groupB.mol(j)]>temp){
                                buff_apos2[m][groupB.mol(j)] = temp;
                                buff_aresid[m][groupB.mol(j)] =groupA.resnum(i);
                                buff_aatmid[m][groupB.mol(j)] =groupA.atom(i);
                                buff_aatom[m][groupB.mol(j)] =groupA.name(i);

                                buff_bmolid[m][groupB.mol(j)] =groupB.mol(j);
                                buff_bresid[m][groupB.mol(j)] =groupB.resnum(j);
                                buff_batmid[m][groupB.mol(j)] =groupB.atom(j);
                                buff_batom[m][groupB.mol(j)] =groupB.name(j);
                              }
                            }
                            if(apos2[m]==0){
                                apos2[m] = temp;
                                aresid[m] =groupA.resnum(i);
                                aatmid[m] =groupA.atom(i);
                                aatom[m] =groupA.name(i);

                                bmolid[m] =groupB.mol(j);
                                bresid[m] =groupB.resnum(j);
                                batmid[m] =groupB.atom(j);
                                batom[m] =groupB.name(j);

                                if(debug>=3)
                                    cout << " apos2 initialized : " << m << " : " << i << "\t" << groupB.mol(j) << " : " << j << "\t" << apos2[m] << endl;
                            }else {
                              if(apos2[m]>temp){
                                if(debug>=3)
                                    cout << " redefine apos2 : " << m << " : " << i << "\t" << groupB.mol(j) << " : " << j << "\t" << temp << endl;
                                apos2[m] = temp;
                                //cout << " test 1 " << aresid[m] << "\t" <<groupA.resnum(i) << endl;
                                aresid[m] =groupA.resnum(i);
                                //cout << " test 2 " << aatom[m] << "\t" <<groupA.name(i) << endl;
                                aatmid[m] =groupA.atom(i);
                                aatom[m] =groupA.name(i);

                                //cout << " test 3 " << bmolid[m] << "\t" <<groupB.mol(j) << endl;
                                bmolid[m] =groupB.mol(j);
                                //cout << " test 4 " << aresid[m] << "\t" <<groupB.resnum(j) << endl;
                                bresid[m] =groupB.resnum(j);
                                //cout << " test 5 " << aatom[m] << "\t" <<groupB.resnum(j) << endl;
                                batmid[m] =groupB.atom(j);
                                batom[m] =groupB.name(j);
                                //cout << " test 6 " << endl;

                                if(apos2[m]<dc ){
                                    pass=m;
                                    break;
                              }
                            }
                        }
                    }else{
                        if(debug>=3)
                            cout << "# this is an H atom, thus skipped. " << endl;
                    }
            }
          }
	}
        for(unsigned int i=0;i<buff_apos2.size();++i){
            for(unsigned int j=0;j<buff_apos2[i].size();++j)
                  if(buff_apos2[i][j]!=0)
                      cout << "# " << setw(8) << numFrames
                                            << setw(5) << i+1
                                            << ":res(" << setw(5) << buff_aresid[i][j]+1
                                            << ":"     << setw(5) << buff_aatom[i][j]
                                            << ") "    << setw(5) << buff_bmolid[i][j]+1
                                            << ":res(" << setw(5) << buff_bresid[i][j]+1
                                            << ":"     << setw(5) << buff_batom[i][j]
                                            << ") # "  << setw(6) << buff_aatmid[i][j]+1
                                            << setw(8) << buff_batmid[i][j]+1
                                            << setw(12) << sqrt(buff_apos2[i][j]) << endl;
        }

        cout << endl << "# Summary" << endl <<  "# =======" << endl
            << setw(10) << "#    Frame" << setw(8) << "molA"
            << setw(8) << "resA" << setw(6) << "atomA"
            << setw(8) << "molB" << setw(8) << "resB" << setw(6) << "atomB"
            << setw(10) << "atomA" << setw(8) << "atomB"
            << setw(12) << "dist [nm]" << endl;
        for(unsigned int i=0;i<apos2.size();++i){
            if(apos2[i]!=0.0)
                cout << "# " << setw(8) << numFrames
                        << setw(5) << i+1
                        << ":res(" << setw(5) << aresid[i]+1
                        << ":"     << setw(5) << aatom[i]
                        << ") "    << setw(5) << bmolid[i]+1
                        << ":res(" << setw(5) << bresid[i]+1
                        << ":"     << setw(5) << batom[i]
                        << ") # "  << setw(6) << aatmid[i]+1
                        <<  setw(8) << batmid[i]+1
                        << setw(12) << sqrt(apos2[i]) << endl;
        }

        numFrames++;
        if(debug){
            cout << endl << "# time used (s) : " << (clock()-clo)/(CLOCKS_PER_SEC) << endl;
        }
        }
      }
      ic.close();
    } //end loop over trajectory
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

