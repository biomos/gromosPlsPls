
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
#include "../src/utils/ExpressionParser.h"

#include "message.h"
#include "block.h"
#include "jobinfo.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

void parse_update
(
 std::vector<std::string> & buffer,
 Jobinfo & jobinfo,
 ExpressionParser<double> & parser,
 std::map<std::string, std::vector<ExpressionParser<double>::expr_struct> > & expr,
 std::map<std::string, data_type> & type,
 std::map<std::string, double> & parameter
);

std::string substitute_var
(
 std::string s, Jobinfo & jobinfo, int j
 );


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
    ExpressionParser<double> parser;
    // map to store ALL input parameters for
    // ALL jobs in the format
    // jobid.block.var
    std::map<std::string, double> parameter;
    std::map<std::string, std::vector<ExpressionParser<double>::expr_struct> > expr;
    std::map<std::string, data_type> type;

    // let's see what we have to do
    Jobinfo jobinfo;
    {
      Ginstream gij(args["jobinfo"]);
      vector<string> buffer;
      gij.getblock(buffer);
      jobinfo.read(buffer);
    }
    
    {
      Ginstream gif(args["format"]);
      vector<string> buffer;
      while (!gif.stream().eof()){
	gif.getblock(buffer);
	
	if(buffer.size()){
	  std::cerr << "reading block " << buffer[0] << std::endl;

	  // checking for special blocks
	  if (buffer[0] == "UPDATE"){
	    parse_update(buffer, jobinfo, parser, expr, type, parameter);
	  }
	  else{
	    Block b;
	    b.init(buffer);
	    input_block.push_back(b);
	  }
	  
	}
      }
    }
    
    // read topology into block
    // read coordinates
    // read special files

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

    
    // check the original input
    {
      Message message;
      
      for(unsigned int i=0; i<input_block.size(); ++i){
	input_block[i].check(message, input_block);
      }
      
      message.display();
      
    }

    // extract once without updating!
    // these will be the values used if
    // PREVID (RUNAFTER) == -1 for a job
    for(unsigned int i=0; i<input_block.size(); ++i){
      input_block[i].extract_var(-1, parameter);
    }
        
    std::cerr << "================================================================================\n"
	      << "  job preparation\n"
	      << "================================================================================\n"
	      << std::endl;
    
    std::map<int, std::map<std::string, int> >::iterator
      j_it = jobinfo.int_data.begin(),
      j_to = jobinfo.int_data.end();

    for( ; j_it != j_to; ++j_it){
      
      const int j = j_it->first;
      if (j != j_it->second["JOBID"])
	throw gromos::Exception("genjob", "JOBID does not match map key");
      
      std::cerr << "JOB " << jobinfo.int_data[j]["JOBID"] << std::endl;

      // update
      for(unsigned int i=0; i<input_block.size(); ++i){
	std::cerr << "update " << input_block[i].name << std::endl;
	
	input_block[i].update(jobinfo.int_data[j],
			      jobinfo.real_data[j],
			      jobinfo.word_data[j]);
	
	input_block[i].extract_var(jobinfo.int_data[j]["RUNAFTER"], parameter);
      }
    } // loop over jobs
    
    // do the expressions from the UPDATE block
    parser.calculate(expr, parameter);

    // back-substitute the calculated values into the jobinfo maps
    jobinfo.back_substitute(parameter, type);

    std::cerr << "================================================================================\n"
	      << "  job check and write\n"
	      << "================================================================================\n"
	      << std::endl;

    j_it = jobinfo.int_data.begin();
    for(; j_it != j_to; ++j_it){

      const int j = j_it->first;
      if (j != j_it->second["JOBID"])
	throw gromos::Exception("genjob", "JOBID does not match map key");

      std::cerr << "JOB " << jobinfo.int_data[j]["JOBID"] << std::endl;
      
      // check
      Message message;
      
      for(unsigned int i=0; i<input_block.size(); ++i){
	input_block[i].check(message, input_block);
      }
      
      
      // update again 
      // (this time using also the values calculated in the UPDATE block)
      for(unsigned int i=0; i<input_block.size(); ++i){
	std::cerr << "update " << input_block[i].name << std::endl;
	
	input_block[i].update(jobinfo.int_data[j],
			      jobinfo.real_data[j],
			      jobinfo.word_data[j]);
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
  catch (std::runtime_error e){
    std::cerr << "RUNTIME ERROR:\n\t";
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

void parse_update(std::vector<std::string> & buffer,
		  Jobinfo & jobinfo,
		  ExpressionParser<double> & parser,
		  std::map<std::string, std::vector<ExpressionParser<double>::expr_struct> > & expr,
		  std::map<std::string, data_type> & type,
		  std::map<std::string, double> & parameter)
{

  for(unsigned int i=1; i<buffer.size()-1; ++i){
    
    // each line should have the following format
    // block.name:TYPE = expression
    // everything in between '$' signs will be replaced by
    // jobid parameters and the :TYPE will be stripped and
    // stored seperately.
    
    std::cerr << "trying to add expression: " << buffer[i] << std::endl;

    std::string::size_type it = buffer[i].find('=');
    std::string::size_type col = buffer[i].find(':');
    
    if (it == std::string::npos || col == std::string::npos ||
	col > it)
      throw gromos::Exception("UPDATE", "wrong 'name' format in UPDATE block");
    
    // now the first part is the variable name, the second the type and
    // the rest the expression
    std::string name;
    {
      std::istringstream is(buffer[i].substr(0, col));
      is >> name;
    }
    
    std::string stype;
    {
      std::istringstream is(buffer[i].substr(col+1, it-col));
      is >> stype;
    }

    // and the whole rest
    std::string e =buffer[i].substr(it+1, std::string::npos);
    std::cerr << "\texpression " << name << " = " << e << std::endl;
    
    
    std::map<int, std::map<std::string, int> >::iterator
      j_it = jobinfo.int_data.begin(),
      j_to = jobinfo.int_data.end();

    for( ; j_it != j_to; ++j_it){

      const int j = j_it->first;
      
      std::string job_name = substitute_var(name, jobinfo, j);
      std::string job_e = substitute_var(e, jobinfo, j);
      
      std::cerr << "\tadding expression " << job_name << " = " << job_e << std::endl;

      std::vector<ExpressionParser<double>::expr_struct> single_expr;
      parser.parse_expression(job_e, parameter, single_expr);
      
      expr[job_name] = single_expr;
      
      std::cerr << "\tparse_update : type == " << stype << std::endl;

      if(stype == "int" || stype == "INT")
	type[job_name] = int_type;
      else if (stype == "real" || stype == "REAL")
	type[job_name] = real_type;
      else if (stype == "word" || stype == "WORD")
	type[job_name] = word_type;
      
      
    }
  }
  
}

std::string substitute_var(std::string s, Jobinfo & jobinfo, int j)
{
  std::string::size_type it = s.find('$');
  if (it == std::string::npos) return s;
  
  std::string::size_type it2 = s.find('$', it+1);
  if (it2 == std::string::npos)
    throw gromos::Exception("substitute_var", "closing '$' not found in '" + s + "'");
  
  std::ostringstream os;
  
  os << s.substr(0, it);

  std::string name = s.substr(it+1, it2 - it - 1);
  std::cerr << "substitute var: " << name << std::endl;

  // special case
  if (name == "PREVID" || name == "RUNAFTER"){
    if (jobinfo.int_data[j]["RUNAFTER"] == -1){
      // insert nothing! and replace the dot
      ++it2;
    }
    else{
      os << "ID" << jobinfo.int_data[j]["RUNAFTER"];
    }
  }
  else if (name == "JOBID" || name == "ID"){
    os << "ID" << jobinfo.int_data[j]["JOBID"];
  }
  else if (jobinfo.int_data[j].count(name))
    os << jobinfo.int_data[j][name];
  else if (jobinfo.real_data[j].count(name))
    os << jobinfo.real_data[j][name];
  else if (jobinfo.word_data[j].count(name))
    os << jobinfo.word_data[j][name];
  else
    throw gromos::Exception("substitute_var", "could not substitute '" + name + "'");
  
  os << s.substr(it2+1, std::string::npos);
  
  return substitute_var(os.str(), jobinfo, j);

}
