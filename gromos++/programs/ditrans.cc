//time series tser

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
// i will use properties
#include "../src/utils/PropertyContainer.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

double nearest_minimum(double phi, double cosdelta, int multiplicity);


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "verbose", "prop", "traj", "strict"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@prop   <property specifier>\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@time   <T and dtn";
  usage += "\t[@strict]\n";
  usage += "\t[@verbose]\n";
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    
    // check if verbose output required
    bool verb = false;
    if (args.count("verbose") >= 0)
      verb = true;
    
    // check if strict checking requested
    bool strict = false;
    if (args.count("strict") >= 0)
      strict = true;

    //  read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    
    System sys(it.system());
    GromosForceField gff(it.forceField());
        
    //   get simulation time
    double time=0, dt=1; 
    {
      Arguments::const_iterator iter=args.lower_bound("time");
      if(iter!=args.upper_bound("time")){
        time=atof(iter->second.c_str());
        ++iter;
      }
      if(iter!=args.upper_bound("time"))
        dt=atof(iter->second.c_str());
    }
  

    // read in a property
    
    // it's nice to actually store the properties read in somewhere
    // let's use a PropertyContainer, as it can handle some of the
    // work later on
    // the PropertyContainer should know about the system, so it can
    // check the input whether ie the atoms exist in the topology
    PropertyContainer props(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("prop");
      Arguments::const_iterator to=args.upper_bound("prop");
      // we read in all properties specified by the user
      // all should be dihedrals!!! (check that???)
      if(iter==to)
	throw Arguments::Exception("no property given");
      for(int i=0; iter!=to; iter++, ++i)
	{
	  string spec=iter->second.c_str();
	  // and that's how easy it is to add a standard property
	  // like distance, angle, torsional angle
	  props.addSpecifier(spec);
	  if (props[i]->type() != "Torsion")
	    throw 
	      Arguments::Exception("Only dihedrals (torsion) properties allowed");
	}    
    }

    // loop over all properties and store the dihedral types
    vector<int> dihedral_types;
    vector<int> number_transitions;
    
    vector<double> old_min;
    
    for(unsigned int i=0; i< props.size(); i++){
      int t=props[i]->getTopologyType(sys);
      if(t<0)
	throw gromos::Exception("ditrans", "Property "+props[i]->toTitle()+
				" not found");
      dihedral_types.push_back(t);
      cout << "# property: " << props[i]->toTitle() << " type " << t << endl;
      old_min.push_back(0.0);
      number_transitions.push_back(0);
      
    }
    //test the nearest minumum function
    //  for(int i=0; i<360; i++){
    // double a=i/180.0*M_PI;
    //int m=6;
    //double cosdelta=-1;
      
    //cout << i 
    ///   << "\t" << 1+cosdelta*cos(m*a) 
    //   << "\t" << 180.0/M_PI*nearest_minimum(a, cosdelta, m) << endl;
    //}
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // define input coordinate
    InG96 ic;
    bool lfirst=true;
    
    // loop over all trajectories
    for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys;
	(*pbc.*gathmethod)();
       
	// calculate the props
	// this is now the place, where a property-container is very handy
	// it knows that we want to loop over all properties and calculate
	// their 'value'. This works best if the property can be represented
	// by one number (which is often the case)
	props.calc();
	if(lfirst){
	  lfirst=false;
	  for(unsigned int i=0; i<props.size(); i++){
	    old_min[i]=nearest_minimum(props[i]->getValue(), 
				      gff.dihedralType(dihedral_types[i]).pd(),
				      gff.dihedralType(dihedral_types[i]).np());
	  }
	}
	else{
	  if (!strict){
	    for(unsigned int i=0; i<props.size(); i++){
	      double min=nearest_minimum(props[i]->getValue(),
					 gff.dihedralType(dihedral_types[i]).pd(),
					 gff.dihedralType(dihedral_types[i]).np());
	      cout << time << "\t" << props[i]->getValue() << "\t" << min << "\t" << old_min[i] << endl;
	      
	      if(min!=old_min[i]){
		if(verb)
		  cout << "# at time " << time << ": transition of " 
		       << props[i]->toTitle() << " from " << old_min[i] 
		       << " to " << min << endl;

		number_transitions[i]++;
		old_min[i]=min;

	      }
	    }
	  }
	  else{
	    for(unsigned int i=0; i<props.size(); i++){

	      double delta_phi = 360.0 / gff.dihedralType(dihedral_types[i]).np();
	      if (abs(old_min[i] - props[i]->getValue()) > delta_phi){
		double min = 
		  nearest_minimum(props[i]->getValue(),
				  gff.dihedralType(dihedral_types[i]).pd(),
				  gff.dihedralType(dihedral_types[i]).np());

		if(verb)
		  cout << "# at time " << time << ": transition of " 
		       << props[i]->toTitle() << " from " << old_min[i] 
		       << " to " << min << endl;
		
		number_transitions[i]++;
		old_min[i]=min;

	      }
	    }	    
	  }
	}
	
	time+=dt;
	
      }
      
    
    
      ic.close();
    }
    for(unsigned int i=0; i< props.size(); i++){
      cout << props[i]->toTitle() << "\t" << number_transitions[i] << endl;
    }
    
  }
  catch (const gromos::Exception &e){
    cerr << "EXCEPTION:\t";
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


double nearest_minimum(double phi, double cosdelta, int multiplicity)
{
  double a_minimum = 180.0 * (3.0-cosdelta) / (2.0*multiplicity);
  double delta_phi = 360.0 / multiplicity;
  //  double minimum_within_pi=a_minimum - int(rint((a_minimum - phi)/delta_phi))*2*M_PI;
  double nearest_min=a_minimum
    - int(rint((a_minimum - phi)/delta_phi))*delta_phi;
  // not so nice to get it down to zero again fi cosdelta = -1 and nearest_min
  // is almost 360.0
  if(nearest_min >= 360.0-0.001) nearest_min-=360.0;
  
  return nearest_min;
}

