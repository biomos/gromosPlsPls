// args_GatherParser.cc

#include <cassert>

#include "GatherParser.h"
#include "Arguments.h"
#include "../gcore/System.h"
#include "../bound/Boundary.h"
#include "../utils/AtomSpecifier.h"
#include "BoundaryParser.h"

using namespace args;
using namespace bound;
using namespace utils;
using namespace gcore;


bound::Boundary::MemPtr GatherParser::parse(gcore::System &sys,gcore::System &refSys,
        const Arguments &gathargs,
        const std::string &str){

    string usage = "\n\nAvailable gathering methods\n";
    usage += "\n\t1 (or glist) : based on a list (default). If no list is given, a warning message will be shown";
    usage += "\n\t2 (or gtime) : based on previous frame (time)";
    usage += "\n\t3 (or gref)  : based on a reference structure (refg required)";
    usage += "\n\t4 (or gltime): the 1st frame based on a list, then previous frame (list+time)";
    usage += "\n\t5 (or grtime): the 1st frame based on a ref, then previous frame (ref+time, refg required)";
    usage += "\n\t6 (or gbond) : based on bond connectivity (useful for a single molecule)";
    usage += "\n\t0 (or nog)   : no gathering\n";

    usage += "\nTwo sub-options are available:";
    usage += "\n\t list   [atom_list (optional for methods 1 and 4)]";
    usage += "\n\t refg   <reference molecule for gathering (required for methods 3 and 5)>\n";

    Boundary::MemPtr gathmethod;

  try {
    Arguments::const_iterator it = gathargs.lower_bound(str);
    Arguments::const_iterator to = gathargs.upper_bound(str);
    if (it == gathargs.upper_bound(str))
      throw Arguments::Exception("###### GATHER WARNING ######\n" + usage);
    ++it;
    // define bool for checking
    int l=0; // whether a list is available
    int r=0; // whether a ref is available

    if (it == gathargs.upper_bound(str)) {
      gathmethod = &Boundary::gatherlist;
      std::cout << "###### GATHER WARNING ######\n"
              << "# Gathering : You have requested to gather the system based on \n"
              << "# an atom list, while you didn't define such a list, therefore \n"
              << "# the gathering will be done according to the 1st atom of the previous molecule\n" << endl;
    } else {
      std::string gather = it->second;
      cout << "# gather option : " << gather << endl;

      if (gather == "nog" || gather == "0") {
        gathmethod = &Boundary::nogather;
      } else if (gather == "1" || gather == "glist") {
        gathmethod = &Boundary::gatherlist;
      } else if (gather == "2" || gather == "gtime") {
        gathmethod = &Boundary::gathertime;
      } else if (gather == "3" || gather == "gref") {
        gathmethod = &Boundary::gatherref;
      } else if (gather == "4" || gather == "gltime") {
        gathmethod = &Boundary::gatherltime;
      } else if (gather == "5" || gather == "grtime") {
        gathmethod = &Boundary::gatherrtime;
      } else if (gather == "6" || gather == "gbond") {
        gathmethod = &Boundary::gatherbond;
      } else {
        throw gromos::Exception("Gather", gather +
                " unknown. " + usage );
      }

      ++it;

      //gcore::System sys(*pbc.sys());
      //System refSys(*pbc.refSys());

      if (it == gathargs.upper_bound(str)) {
          if(gather=="1" || gather == "4"){
              cout << "# Gathering : You have requested to gather the system based on " << endl
                    << "# an atom list, while you didn't define such a list, therefore "<< endl
                    << "# the gathering will be done according to the 1st atom of the previous molecule" << endl;
                }
          if(gather=="3" || gather == "5"){
              throw gromos::Exception("Gathering : ",gather +
                      " No reference structure given!");
              }
      } else {
          AtomSpecifier gathlist(sys);

//          string gathopt = it->second;
          //for(;it!=to;it++){
          while(it!=to){
              string gathopt = it->second;
              if(gathopt=="list"){
                  it++;
                  while(gathopt!="refg"&&it!=to){
                      string gathopt=it->second.c_str();
                      gathlist.addSpecifierStrict(gathopt);
                      it++;
                      //cout << "# gather list : " << gathopt << endl;
                  }
                  //the block of primlist
                  if(gather=="1" || gather=="4")
                  for(int j=0;j<gathlist.size()/2;++j){
                        int i=2*j;
                        sys.primlist[gathlist.mol(i)][0]=gathlist.atom(i);
                        sys.primlist[gathlist.mol(i)][1]=gathlist.mol(i+1);
                        sys.primlist[gathlist.mol(i)][2]=gathlist.atom(i+1);

                        refSys.primlist[gathlist.mol(i)][0]=gathlist.atom(i);
                        refSys.primlist[gathlist.mol(i)][1]=gathlist.mol(i+1);
                        refSys.primlist[gathlist.mol(i)][2]=gathlist.atom(i+1);

                        cout << "# updated prim : mol " << gathlist.mol(i) << " atom " << gathlist.atom(i)
                                << "\t# refe : mol " << sys.primlist[gathlist.mol(i)][1]
                                << " atom " << sys.primlist[gathlist.mol(i)][2] << endl;
                    }
                  l=1;
              } else if(gathopt=="refg") {
                  it++;
                  if(it==to){
                      throw gromos::Exception("Gathering : ",gather +
                      " No reference structure given!");
                  }

                  if(gather=="3" || gather=="5"){
                  string reffile = it->second.c_str();

                  Boundary *pbc = BoundaryParser::boundary(sys, gathargs);
                  pbc->setReferenceFrame(reffile);
                  //cout << "# gather ref : " << reffile << endl;
                  }

                  // in case refg comes before list
                  it++;
                  r=1;

              } else {
                  throw gromos::Exception("Unknown gather option ",gathopt +
                          " \nKnown gather options are : list + refg!");
              }
          }
      }
      if(l==0 && (gather=="1" || gather=="4")){
          cout << "\n# Gathering : You have requested to gather the system based on " << endl
               << "\n# an atom list, while you didn't define such a list, therefore "<< endl
               << "\n# the gathering will be done according to the 1st atom of the previous molecule" << endl;
      }
      if(r==0 && (gather=="3" || gather=="5")){
          throw gromos::Exception("Gathering : ",gather +
                      " No reference structure given!");
      }
    }
  } catch (Arguments::Exception &e) {
    //gathmethod = &Boundary::coggather;
      gathmethod = &Boundary::gatherlist;
  }


  return gathmethod;
}

