// checkbox.cc
//
// Program to calculate, for every frame in a set of trajectories, the
// minimum distance between two periodic copies of a (gathered) molecule
//
// notes: only works for truncated octahedral boundary conditions "-pbc t"
//        do not know if gather can gather multiple molecules correctly
//
// adapted to pbc r --mika
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gmath/Vec.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/TrajArray.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"


#include <vector>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace fit;
using namespace utils;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "mol", "pbc"};
  int nknowns = 4;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@mol <molecules to include in calculation (default: all)>\n";
  usage += "\t@traj <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    
    // Systems for calculation
    System sys(it.system());

    // Periodic images of system
    const int max_num_of_images=14;
    System *compare_periodic[max_num_of_images];
    for (int ii=0; ii<max_num_of_images; ii++) {
      compare_periodic[ii]= new System(sys);
    }
    gmath::Vec displace[max_num_of_images];

    // Coordinates input file
    InG96 ic;

    // Coordinates temporary storage
    TrajArray ta(sys);

    //get the molecules to be included
    int mol = sys.numMolecules();
       {
   Arguments::const_iterator iter=args.lower_bound("mol");
   if(iter!=args.upper_bound("mol")){
     mol=atoi(iter->second.c_str());
   }
    }
       if (mol > sys.numMolecules()) {throw gromos::Exception(
    "checkbox", " cannot go beyond the total number of molecules in topology!\n");
    }

    // Parse boundary conditions for sys
       int num_of_images=0;
       {
	 Arguments::const_iterator iter=args.lower_bound("pbc");
	 char b=iter->second.c_str()[0];
	 switch(b){
	 case 't': num_of_images=14;
	 case 'r': num_of_images=6;
	 }
       }
       
      Boundary *pbc = BoundaryParser::boundary(sys, args);

   // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);


    // Distance calculation variables
    gmath::Vec distvec;
    double overall_min_dist2=9.9e12;
    double mindist2=9.9e12;

    int frame=0;
    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);

      // loop over all frames
      while(!ic.eof()){
	ic >> sys; frame++;

     (*pbc.*gathmethod)();
        ta.store(sys,0); 
        for (int ii=0; ii<num_of_images; ii++) {
          ta.extract(*compare_periodic[ii],0);
          for (int jj=0; jj<3; jj++) displace[ii][jj]=sys.box()[jj];
        }
        if (num_of_images >= 6) {
          displace[ 0][0] *= 1.0; displace[ 0][1] *= 0.0; displace[ 0][2] *= 0.0;
          displace[ 1][0] *=-1.0; displace[ 1][1] *= 0.0; displace[ 1][2] *= 0.0;
          displace[ 2][0] *= 0.0; displace[ 2][1] *= 1.0; displace[ 2][2] *= 0.0;
          displace[ 3][0] *= 0.0; displace[ 3][1] *=-1.0; displace[ 3][2] *= 0.0;
          displace[ 4][0] *= 0.0; displace[ 4][1] *= 0.0; displace[ 4][2] *= 1.0;
          displace[ 5][0] *= 0.0; displace[ 5][1] *= 0.0; displace[ 5][2] *=-1.0;
        }
        if (num_of_images == 14) {
          displace[ 6][0] *= 0.5; displace[ 6][1] *= 0.5; displace[ 6][2] *= 0.5;
          displace[ 7][0] *= 0.5; displace[ 7][1] *= 0.5; displace[ 7][2] *=-0.5;
          displace[ 8][0] *= 0.5; displace[ 8][1] *=-0.5; displace[ 8][2] *= 0.5;
          displace[ 9][0] *= 0.5; displace[ 9][1] *=-0.5; displace[ 9][2] *=-0.5;
          displace[10][0] *=-0.5; displace[10][1] *= 0.5; displace[10][2] *= 0.5;
          displace[11][0] *=-0.5; displace[11][1] *= 0.5; displace[11][2] *=-0.5;
          displace[12][0] *=-0.5; displace[12][1] *=-0.5; displace[12][2] *= 0.5;
          displace[13][0] *=-0.5; displace[13][1] *=-0.5; displace[13][2] *=-0.5;
        }
        for (int ii=0; ii<num_of_images; ii++) {
          PositionUtils::translate(compare_periodic[ii],displace[ii]);
        }

        mindist2=9.9e12;

        for(int m1=0;m1<mol;++m1) {
          for(int m2=0;m2<mol;++m2) {
            for (int i1=0;i1 < sys.mol(m1).numAtoms(); i1++) {
              for (int i2=0;i2 < sys.mol(m2).numAtoms(); i2++) {
                for (int image=0; image<num_of_images; image++) {
                  distvec= sys.mol(m1).pos(i1) 
                          -  compare_periodic[image]->mol(m2).pos(i2);
                  if (distvec.abs2() < mindist2) {
                    mindist2=distvec.abs2();
                  }
                }
              }
            }
          }
        }
        cout << setw(11) << frame;
        cout.precision(7);
        cout << setw(15) << sqrt(mindist2) << endl;

        if (mindist2 < overall_min_dist2) {
          overall_min_dist2=mindist2;
        }

      }
      ic.close();
    }
    cout << "OVERALL MIN";
    cout.precision(7);
    cout << setw(15) << sqrt(overall_min_dist2) << endl;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

