
#include <cassert>

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include "../gmath/Vec.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"

#include "PropertyContainer.h"

using namespace gcore;
using namespace std;
using namespace utils;

namespace utils
{
  typedef vector<Property *>::value_type argType;
  
  //---PropertyContainer Class------------------------------------

  PropertyContainer::PropertyContainer() :
    vector<Property *>()
  {
    d_sys = NULL;
    d_distribution = NULL;
  }
  
  PropertyContainer::PropertyContainer(gcore::System &sys) :
    vector<Property *>()
  {
    d_sys = &sys;
    d_distribution = NULL;
  }
  
  PropertyContainer::~PropertyContainer()
  {
  }
  
  int PropertyContainer::addSpecifier(std::string s)
  {
    parse(s);
    return size();
  }
  
  void PropertyContainer::parse(std::string s)
  {
    std::string::size_type iterator = 0, last = 0;
    std::string arguments[Property::MAXARGUMENTS];
    std::string type = "";
    int count = 0;
        
    // extract type    
    iterator = s.find('%', iterator);
    if (iterator == std::string::npos)
      throw PropertyContainer::Exception(
            " invalid property-specifier.\nSyntax: <type>\%<atomspecifier>[\%...]\n");
    
    type = s.substr(0, iterator); // allows for arbitrary string length types

    // now separate the string on % into argument tokens
    for(; count < Property::MAXARGUMENTS; ++count)
      {
	last = iterator;
	iterator = s.find('%', iterator+1);
	if (iterator == std::string::npos)
	  {
	    // the last argument
	    arguments[count] = s.substr(last+1, s.length());
	    // just to the end of the string
	    ++count;
	    break;
	  }
	arguments[count] = s.substr(last+1, iterator-last-1);
      }

    // check if the property should be added for all molecules
    // if the fisrt argument (after the type) is not of standard
    // format, it must not begin with a 'a:'!!!
    // otherwise, the a will get substituted with all molecule numbers
    // one at a time...

    if (count >=1 && (arguments[0]).substr(0,2) == "a:")
      {
	// add this property for all molecules
	// issue a warning, if more molecules are specified
	if ((arguments[0]).find(':',2) != std::string::npos)
	  cerr << "# WARNING: adding a _inter_ molecular property"
	       << " for all molecules!\n";
	// loop over d_sys->numMolecules()
	// exchange substring 'a:' with the numbers of the molecules
	// insert into container
	// cerr << "# adding all molecules\n";
	char b[100];       
	std::string rest = arguments[0].substr(1,arguments[0].length()-1);
	for(int i=0; i<d_sys->numMolecules(); i++)
	  {
	    // g96 arrays start with 1
	    sprintf(b, "%d", i+1);
	    std::string s = b;	    
	    arguments[0] = s + rest;
	    insert(end(), argType(createProperty(type, count, arguments)));
	  }
	return;
	
      }
    // normal insert of one specified property (not all molecules)
    insert(end(), argType(createProperty(type, count, arguments)));
  }

  Property* PropertyContainer::createProperty(std::string type, int count, 
				     std::string arguments[])
  {
    // this method has to know all existing property types
    // when adding user properties, a child class of PropertyContainer
    // has to be added which extends this method for the new properties

    if (type == "d")
      { 
	DistanceProperty *p = new DistanceProperty(*d_sys);
	p->parse(count, arguments);
	return p;
      }
    
    if (type == "a")
      {
	AngleProperty *p = new AngleProperty(*d_sys);
	p->parse(count, arguments);
	return p;
      }
    
    if (type == "t"||type == "i")
      {
	TorsionProperty *p = new TorsionProperty(*d_sys);
	p->parse(count, arguments);
	return p;
      }

    if (type == "o")
      {
	OrderProperty *p = new OrderProperty(*d_sys);
	p->parse(count, arguments);
	return p;
      }

    // throw exception type error
    // either do the user properties first or catch the exception
    throw PropertyContainer::Exception(" type unknown\n");

    return NULL;
  }

  std::string PropertyContainer::toTitle()
  {
    std::string s = "";
    // stop output after 10 properties
    iterator it = begin();
    for(int c = 0; it != end() && c<10; it++, c++)
      s += (*it)->toTitle() + "\t\t";
    if (it != end())
      cout << "...";
    return s;
  }
  
  std::string PropertyContainer::toString()
  {
    std::string s = "";
    for(iterator it = begin(); it != end(); it++)
      s += (*it)->toString() + "\t\t";
    return s;
  }
  
  void PropertyContainer::calc()
  {
    // distribution can only be used, if property returns one double
    if (d_distribution)
      for(iterator it = begin(); it != end(); ++it)
	d_distribution->add((*it)->calc());
    else
      for(iterator it = begin(); it != end(); ++it)
	(*it)->calc();
  }

  std::string PropertyContainer::averageOverRun()
  {
    std::string s = "";
    for(iterator it = begin(); it != end(); ++it)
      s += (*it)->average() + "\t\t";
    return s;
  }
  
  std::string PropertyContainer::averageOverProperties()
  {
    // gives <average> <rmsd> <rmsd from zvalue> <lowes value> <highest value>
    // averaged over all properties !!!
    // this is intended for 'all molecule' properties
    // it should be called after every calculation
    double av;
    double ub, lb;
    double rmsd, zrmsd;

    averageOverProperties(av, rmsd, zrmsd, lb, ub);
    
    ostringstream os;
    os << av << "\t" << rmsd << "\t"
			  << zrmsd << "\t" << lb << "\t"
			  << ub;
    return os.str();

    // sprintf(b, "%f\t\t%f\t\t%f\t\t%f\t\t%f", av, rmsd, zrmsd, lb, ub);
    // std::string s = b;
    // return s;
  }

  void PropertyContainer::averageOverProperties(double &av, double &rmsd, 
						double &zrmsd, double &lb, 
						double &ub)
  {
    // gives <average> <rmsd> <rmsd from zvalue> <lowes value> <highest value>
    // averaged over all properties !!!
    // this is intended for 'all molecule' properties
    // it should be called after every calculation
    double v;
    av = 0.0;
    ub = -MAXFLOAT;
    lb = MAXFLOAT;
    zrmsd = 0.0;
    
    for(iterator it = begin(); it != end(); ++it)
      {
	v = (*it)->getValue();
	av += v;
	if (v > ub) ub = v;
	if (v < lb) lb = v;
	zrmsd += pow(v-(*it)->getZValue(), 2);
      }
    av /= size();
    zrmsd /= size();
    zrmsd = sqrt(zrmsd);
    
    rmsd = 0;
    for(iterator it = begin(); it != end(); ++it)
      rmsd += pow((*it)->getValue() - av, 2);
    rmsd /= size();
    rmsd = sqrt(rmsd);
  }
  
  std::ostream &operator<<(std::ostream &os, PropertyContainer &ps)
  {
    for(PropertyContainer::iterator it = ps.begin(); it != ps.end(); it++)
      os << (**it) << "\t\t";
    os << endl;

    return os;
  }
  
}


  











