
#include <cassert>
#include <sstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gio/Ginstream.h"

#include "block.h"
#include "jobinfo.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv){
  
  char *knowns[] = { "format", "input", "jobinfo"
  };
  
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@format    <input file specification>";
  usage += "\n\t@input     <input file>";
  usage += "\n\t@jobinfo   <jobinfo file>";
  usage += "\n\n";
    
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    
    std::vector<Block> input_block;

    {
      Ginstream gif(args["format"]);
      vector<string> buffer;
      while (!gif.stream().eof()){
	gif.getblock(buffer);
	
	if(buffer.size()){
	  std::cerr << "reading block " << buffer[0] << std::endl;
	  Block b;
	  b.init(buffer);
	  input_block.push_back(b);
	}
      }
    }
    
    {
      Ginstream gii(args["input"]);
      vector<string> buffer;
      while (!gii.stream().eof()){
	gii.getblock(buffer);
	if (buffer.size()){
	  for(unsigned int i=0; i<input_block.size(); ++i){
	    input_block[i] << buffer;
	  }
	}
      }
    }
    
    Jobinfo jobinfo;
    {
      Ginstream gij(args["jobinfo"]);
      vector<string> buffer;
      gij.getblock(buffer);
      jobinfo.read(buffer);
    }
    
    for(int j=0; j<jobinfo.size(); ++j){
      
      std::cerr << "JOB " << jobinfo.int_data[j]["JOBID"] << std::endl;

      // update
      for(unsigned int i=0; i<input_block.size(); ++i){
	std::cerr << "update " << input_block[i].name << std::endl;
	
	input_block[i].update(jobinfo.int_data[j],
			      jobinfo.real_data[j],
			      jobinfo.word_data[j]);
      }
      
      // check
      Message message;
      
      for(unsigned int i=0; i<input_block.size(); ++i){
	input_block[i].check(message, input_block);
      }
      
      
      // write the input file
      for(unsigned int i=0; i<input_block.size(); ++i){
	input_block[i].write();
      }
      
      
      message.display();

    } // loop over jobs
    
    
  }
  catch (gromos::Exception e){
    std::cerr << "EXCEPTION:\t";
    std::cerr << e.what() << std::endl;
    return 1;
  }
  
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// helper functions
////////////////////////////////////////////////////////////////////////////////

variable parse_format(std::string w)
{
  std::string::size_type f = w.find(':');
  if (f == std::string::npos)
    throw gromos::Exception("format", "could not parse variable " + w);
  
  std::string w1 = w.substr(0, f);
  std::string w2 = w.substr(f+1, std::string::npos);
	  
  if (w2 == "INT" || w2 == "int")
    return variable(w1, int_type);
  else if (w2 == "REAL" || w2 == "real")
    return variable(w1, real_type);

  throw gromos::Exception("format", "could not parse variable type " + w2 + " of " + w);

}

void parse_data(std::istringstream & is, variable & data,
		std::map<std::string, int> & int_data,
		std::map<std::string, double> & real_data,
		std::map<std::string, std::string> & word_data)
{
  switch(data.type){
    case int_type:
      {
	int d;
	if (!(is >> d))
	  throw gromos::Exception(data.name, "could not read integer for " + data.name);
	int_data[data.name] = d;
	std::cerr << "read int " << d << std::endl;
	break;
      }
    case real_type:
      {
	double d;
	if (!(is >> d))
	  throw gromos::Exception(data.name, "could not read real for " + data.name);
	real_data[data.name] = d;
	std::cerr << "read real " << d << std::endl;
	break;
      }
    case word_type:
      {
	std::string d;
	if (!(is >> d))
	  throw gromos::Exception(data.name, "could not read word for " + data.name);
	word_data[data.name] = d;
	std::cerr << "read word " << d << std::endl;
	break;
      }
    case new_line:
      {
	break;
      }
  }
}
