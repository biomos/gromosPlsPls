// buildbox.cc

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include <vector>
#include <set>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace args;
using namespace std;

class point
{
public:
  point(int i, int j, int k)
  {
    xi=i;
    yi=j;
    zi=k;
  }
  int xi;
  int yi;
  int zi;
  bool const operator<(point const & p)const
  {
    if(xi<p.xi) return true;
    else if(xi==p.xi && yi < p.yi) return true;
    else if(xi==p.xi && yi==p.yi && zi < p.zi) return true;
    return false;
  }
};


int main(int argc, char **argv){

  char *knowns[] = {"topo1", "insx1", "topo2", "insx2", "nsm", "densit", "fraction"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo1 <topology1>\n";
  usage += "\t@insx1 <coordinates for a single molecule of type 1>\n";
  usage += "\n\t@topo2 <topology2>\n";
  usage += "\t@insx2 <coordinates for a single molecule of type 2>\n";
  usage += "\t@nsm <number of molecules per dimension>\n";
  usage += "\t@densit <density of liquid (kg/m^3)>\n";
  usage += "\t@fraction <mole fraction of 1>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    // set some values
    args.check("nsm",1);
    Arguments::const_iterator iter=args.lower_bound("nsm");
    int nsm3=atoi(iter->second.c_str());
    int nsm=nsm3*nsm3*nsm3;

    double densit=atof(args["densit"].c_str());
    double fraction=atof(args["fraction"].c_str());
    int nsm1=int(fraction*nsm);
    int nsm2=nsm-nsm1;
    
    // read topology
    args.check("topo1",1);
    InTopology it1(args["topo1"]);
    System smol1(it1.system());

    
    // read topology
    args.check("topo2",1);
    InTopology it2(args["topo2"]);
    System smol2(it2.system());

    // and calculate some more values
    double weight1=0, weight2=0;
    for(int i=0; i<smol1.numMolecules();i++)
      for(int j=0; j< smol1.mol(i).numAtoms();j++)
        weight1+=smol1.mol(i).topology().atom(j).mass();
    for(int i=0; i<smol2.numMolecules();i++)
      for(int j=0; j< smol2.mol(i).numAtoms(); j++)
	weight2+=smol2.mol(i).topology().atom(j).mass();
    double weight=nsm1*weight1 + nsm2*weight2;
    
    double vtot=(weight*1.66056)/densit;
    double box=pow(vtot,1.0/3.0);
    double box3=box/nsm3;
    Vec box32(box3/2.0, box3/2.0, box3/2.0);
    //cout << "total number of molecules " << nsm << endl;
    //cout << "number of type 1 " << nsm1 << endl;
    //cout << "number of type 2 " << nsm2 << endl;
    
    // read singe atom coordinates...
    InG96 ic;
    args.check("insx1",1);
    ic.open(args["insx1"]);
    ic >> smol1;
    ic.close();
    ic.open(args["insx2"]);
    ic >> smol2;
    ic.close();
        
    Vec rc1=PositionUtils::com(smol1)-box32;
    Vec rc2=PositionUtils::com(smol2)-box32;
    
    PositionUtils::translate(&smol1, -rc1);
    PositionUtils::translate(&smol2, -rc2);

     // new system
    System sys;
    //for(int i=0;i<3;i++){
    //  double *tmp = (double *) &sys.box()[i];
    //  *tmp = box;
    //}
    //initialize random seed
    srand(int(1000*densit));
   
    // set up a grid
    set<point> grid;
    for(int i=0; i<nsm3; i++)
      for(int j=0; j<nsm3; j++)
	for(int k=0; k<nsm3; k++)
	  grid.insert(point(i,j,k));
    
    //first, we do molecule 1
    for(int i=0; i< nsm1; i++){
      //get an random number
      int r=rand();
      long int itry=(r*grid.size())/RAND_MAX;
      set<point>::iterator iter=grid.begin();
      for(int i=0; i<itry; i++) iter++;
      if(iter==grid.end())
	throw gromos::Exception("bin_box", "unlikely but true: "
				"not enough grid points for random number");
    
      //cout << "itry " << itry << "\t" << ipos << "\t" << ix << "\t" << iy << "\t" << iz << endl;

      Vec shift(iter->xi*box3, iter->yi*box3, iter->zi*box3);
      PositionUtils::translate(&smol1, shift);
      for(int k=0; k< smol1.numMolecules(); k++)
	sys.addMolecule(smol1.mol(k));
      PositionUtils::translate(&smol1, -shift);
      grid.erase(iter);
    }
    //then, we do molecule 2
    if(int(grid.size())!=nsm2)
      throw gromos::Exception("bin_box", "the number of grid points left "
			      "after adding the first species is not the "
			      "same as required for the second species");
    // just fill in the remaining grid points
    for(set<point>::iterator iter=grid.begin(), to=grid.end();
	iter!=to; ++iter){
      
      Vec shift(iter->xi*box3, iter->yi*box3, iter->zi*box3);
      PositionUtils::translate(&smol2, shift);
      for(int k=0; k< smol2.numMolecules(); k++)
	sys.addMolecule(smol2.mol(k));
      PositionUtils::translate(&smol2, -shift);
    }
    for(int i=0; i<3; i++)    
      sys.box()[i]=box;
    
    
    // Print the new set to cout
    OutG96S oc;
    ostringstream os;
    os << "Bin_box generated a binary mixture of" << endl;
    os <<  nsm1 << " x "<<args["insx1"] << endl;
    os <<  nsm2 << " x "<<args["insx2"] << endl;
    os << "Density : " << densit << " kg/m^3" << endl;
    
    oc.open(cout);
    oc.writeTitle(string(os.str()));
    oc << sys;
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




