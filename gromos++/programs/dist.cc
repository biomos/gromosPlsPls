//distributions dist

// written by Mika, modified by Markus

// this program takes 'advanced' use of Properties
// if you don't know about them, read first tser.cc
// for standard usage...

#include "../src/args/Arguments.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
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

// a non-standard property is defined here
// standard ones are: distance, angle, torsion-angle
// this one should check, whether the value is in a certain range
// (ie, whether the dihedral angle is gauche or anti)
class MinimumProperty : public Property
{
public:
  // this property carries a distribution of it's values
  // but here we want a distribution over all properties specified
  // by the user (so they should be of similar type - to get a
  // meaningful distribution). So the distribution objet is stored in
  // a PropertyContainer (actually a derived class from PropertyContainer)
  MinimumProperty(gmath::Distribution &dist);
  virtual ~MinimumProperty();
  
  // the parse function get's called to read in user input (@prop ...)
  // it has to be overwritten if any non standard arguments are needed
  void parse(int count, std::string arguments[]);
  
  // methods
  virtual float calc(); // returns the integral
  virtual std::string average(); // might return the average value (dependent on grid size)
  // meaning: all values inside a bin are treated as if they are exactly in the middle,
  // which of course is not true!
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
  // this is the type of the property, 
  // it is mainly needed for diagnostic output
  d_type = "Minimum";
  // specify how many arguments are required
  // this is especially important, if the standard parse() function
  // is used, if it is overwritten anyway, it's your task to ensure this
  REQUIREDARGUMENTS = 2;
  
  // the distribution to add the calculated values to. here: it's in a
  // PropertyContainer derived class (that also holds the properties)
  d_dist = &dist;
  // the range of the 'well' (of the property-values considered 'in')
  d_begin = 0;
  d_end = 0;
}

MinimumProperty::~MinimumProperty()
{
}

void MinimumProperty::parse(int count, std::string arguments[])
{
  // you see, I use REQUIREDARGUMENTS, that's the same as the Property::parse()
  if (count < REQUIREDARGUMENTS)
    throw MinimumProperty::Exception(" need two arguments to define a minimum.\n");
  
  // and just parse the two arguments into the member variales
  // it's still quite easy to get user input...
  if (sscanf(arguments[0].c_str(), "%f", &d_begin) != 1 || sscanf(arguments[1].c_str(), "%f", &d_end) != 1)
    throw MinimumProperty::Exception(" arguments format error.\n");
}

float MinimumProperty::calc()
{
  // go through the distribution, check how many values are
  // in the well
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
  // print out the property
  char b[100];
  sprintf(b, "min(%f-%f)\t%f", d_begin, d_end, d_value);
  std::string s = b;
  return s;
}

std::string MinimumProperty::toTitle()
{
  // print a title line
  char b[100];
  sprintf(b, "min%%(%f-%f)", d_begin, d_end);
  std::string s = b;
  return s;
}

std::string MinimumProperty::average()
{
  // average might or might not make sense...
  // actually I think it does, maybe should implement it
  throw MinimumProperty::Exception(" average: not implemented.\n");
}


// now we have added a new property, so our property-container has to
// know about it to create the appropriate class.
// there are possiblities to use (or implement) dynamical class creation,
// but i rather choose to overwrite the construction method
// this of course can get problematic, if different PropertyContainer exist,
// which all do not know really all the properties, but since here it's
// fairy linear usage, this problem should not occur...
class MinPropertyContainer : public PropertyContainer
{
public:
  // first special point: we need a distribution!
  MinPropertyContainer(gmath::Distribution &dist);
  virtual ~MinPropertyContainer();
  
protected:
  // now this is the method that actually creates the appropriate property-
  // classes, we have to override to create our new property
  virtual Property* createProperty(std::string type, int count, std::string arguments[]);
  
  // this is a second distribution, not the same as
  // the standard property container holds...
  // the imp is really a hack!
  // (this is because we want to use functionality from the standard
  // implementation, but then the standard distribution will also be used.
  // and of course not in the way we like to use it here!)
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
  // so you can see here, that type is not limited to a character
  // it may be longer
  // so here we handle the new property
  if (type == "min")
    {
      MinimumProperty *p = new MinimumProperty(*d_dist);
      p->parse(count, arguments);
      return p;
    }
  // all default ones are handled nicely by the default implementation
  // actually i should deny any other properties to be specified here
  // cause this is a very special container, not used in the normal way
  PropertyContainer::createProperty(type, count, arguments);
  return NULL;
}

// as you can see, it's fairly easy to add a new property


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
  // these are the standard properties we want to calculate
  // at every timestep
  // these will be added to the standard PropertyContainer distribution
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

  // and we are additionaly interested, how many times the property
  // fall into a certain range in the distribution

  // the normal way would be to subclass distribution, but i really wanted
  // to show how to add a new property... ;-)
  // (to use this new distribution, you would of course then also subclass
  // the propertyContainer, but still you would get a cleaner program than this)
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
  Boundary *pbc = BoundaryParser::boundary(sys, args); 
  // parse gather method
  Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

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
      (*pbc.*gathmethod)();
      props.calc();
      cout << props.checkBounds();
      
    }
  }
  ic.close();
  // this already has been the main program!!!

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

  // here the 'integrals' in the specified wells are calculated
  mins.calc();
  cout << "#\t" << mins << endl;
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}





