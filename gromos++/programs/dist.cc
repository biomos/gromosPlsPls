//distributions dist

#include "../src/args/Arguments.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Distribution.h"
#include "../src/utils/PropertyContainer.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

class MinimumProperty : public Property
{
public:
  MinimumProperty(gmath::Distribution &dist);
  virtual ~MinimumProperty();
  
  void parse(int count, std::string arguments[]);
  
  // methods
  virtual float calc(); // returns the integral
  virtual std::string average(); // might return the average value (dependent on grid size)
  virtual std::string toString();
  virtual std::string toTitle();
  
  struct Exception: public gromos::Exception
  {
    Exception(const string &what): gromos::Exception("MinimumProperty", what) {}
  };
  
protected:
  float d_begin, d_end;
  gmath::Distribution *d_dist;
  
};

MinimumProperty::MinimumProperty(gmath::Distribution &dist) :
  Property()
{
  d_type = "Minimum";
  REQUIREDARGUMENTS = 2;
  
  d_dist = &dist;
  d_begin = 0;
  d_end = 0;
}

MinimumProperty::~MinimumProperty()
{
}

void MinimumProperty::parse(int count, std::string arguments[])
{
  if (count < REQUIREDARGUMENTS)
    throw MinimumProperty::Exception(" need two arguments to define a minimum.\n");
  
  if (sscanf(arguments[0].c_str(), "%f", &d_begin) != 1 || sscanf(arguments[1].c_str(), "%f", &d_end) != 1)
    throw MinimumProperty::Exception(" arguments format error.\n");
}

float MinimumProperty::calc()
{
  d_value=0;
  int count;
  for(int i=0; i<d_dist->nSteps(); i++)
    {
      if (d_dist->value(i) >= d_begin)
	{
	  if (d_dist->value(i) > d_end) break;
	  count = (*d_dist)[i];
	  d_value+=count;
	}  
    }
  d_value /= d_dist->nVal();
  return d_value;
}

std::string MinimumProperty::toString()
{
  char b[100];
  sprintf(b, "min(%f-%f)\t%f", d_begin, d_end, d_value);
  std::string s = b;
  return s;
}

std::string MinimumProperty::toTitle()
{
  char b[100];
  sprintf(b, "min%%(%f-%f)", d_begin, d_end);
  std::string s = b;
  return s;
}

std::string MinimumProperty::average()
{
  throw MinimumProperty::Exception(" average: not implemented.\n");
}

class MinPropertyContainer : public PropertyContainer
{
public:
  MinPropertyContainer(gmath::Distribution &dist);
  virtual ~MinPropertyContainer();
  
protected:
  virtual Property* createProperty(std::string type, int count, std::string arguments[]);
  
  // this is a second distribution, not the same as
  // the standard property container holds...
  // the imp is really a hack!
  gmath::Distribution *d_dist;
};

MinPropertyContainer::MinPropertyContainer(gmath::Distribution &dist)
  : PropertyContainer()
{
  d_dist = &dist;
}

MinPropertyContainer::~MinPropertyContainer()
{
}

Property* MinPropertyContainer::createProperty(std::string type, int count, std::string arguments[])
{
  if (type == "min")
    {
      MinimumProperty *p = new MinimumProperty(*d_dist);
      p->parse(count, arguments);
      return p;
    }
  PropertyContainer::createProperty(type, count, arguments);
  return NULL;
}


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "prop", "dist", "min", "traj"};
  int nknowns = 6;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@dist   <lower and upper boundary and number of steps>\n";
  usage += "\t@prop   <propertyspecifier>\n";
  usage += "\t[@min   <propertyspecifier>]\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  //   get distribution parameters
  double begin=0, end=0;
  int nsteps=0; 
  {
    Arguments::const_iterator iter=args.lower_bound("dist");
    if(iter!=args.upper_bound("dist")){
      begin=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("dist")){
      end=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("dist")){
      nsteps=atoi(iter->second.c_str());
    }     
  }
  
  //  read topology
  args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());

  // get properties into PropertySpecifier
  PropertyContainer props(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("prop");
    Arguments::const_iterator to=args.upper_bound("prop");
    for(; iter!=to; iter++)
      {
	string spec=iter->second.c_str();
	props.addSpecifier(spec);
      }    
  }

  // set up distribution arrays
  gmath::Distribution dist(begin, end, nsteps);
  props.addDistribution(dist);

  MinPropertyContainer mins(dist);
  {
    Arguments::const_iterator iter=args.lower_bound("min");
    Arguments::const_iterator to=args.upper_bound("min");
    for(; iter!=to; ++iter)
      {
	string spec=iter->second.c_str();
	mins.addSpecifier(spec);
      }
  }
  
  // Parse boundary conditions
  Boundary *pbc;
  try{
    char b=args["pbc"].c_str()[0];
    switch(b){
      case 't':
        pbc=new TruncOct(&sys);
        break;
      case 'v':
        pbc=new Vacuum(&sys);
        break;
      case 'r':
        pbc=new RectBox(&sys);
        break;
      default:
        throw gromos::Exception("Boundary", args["pbc"] + 
				" unknown. Known boundaries are t, r and v");
	
    }
  }
  catch(Arguments::Exception &e){
    pbc = new Vacuum(&sys);
  }

  // define input coordinate
  InG96 ic;

  // set pi
  // const double pi = 3.1415926535898;
  
  // loop over all trajectories
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

    // open file
    ic.open((iter->second).c_str());
       // ic.select("ALL");
    
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      pbc->gather();
      props.calc();
      cout << props.checkBounds();
      
    }
  }
  ic.close();
  // print out the distribution, calculate the average and rmsd
  cout << "#" << endl;  
  cout << "# number of values calculated: "   
       << props.getDistribution().nVal() << endl;
  cout << "# average value:               "   
       << props.getDistribution().ave() << endl;
  cout << "# RMSD (from distribution):    "   
       << props.getDistribution().rmsd() << endl;

  cout << "# time\t\t" <<  props.toTitle() << endl;

  props.getDistribution().write(cout);

  mins.calc();
  cout << "#\t" << mins << endl;
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}





