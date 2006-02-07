
#ifndef INCLUDED_JOBINFO_H
#define INCLUDED_JOBINFO_H

#include "block.h"

class Jobinfo
{
public:
  int read(std::vector<std::string> & input)
  {
    if (input.size() < 3)
      throw gromos::Exception("jobinfo", "not enough jobinfo data");
    
    // read data format
    std::istringstream is(input[1]);
    
    while (!is.eof()){
      std::string w;
      is >> w;
      
      if (w == "JOBID" || w == "jobid" || w == "job_id"){
	data.push_back(variable("JOBID", int_type));
      }
      else if (w == "dir" || w == "DIR"){
	data.push_back(variable("DIR", word_type));
      }
      else if (w == "runafter" || w == "RUNAFTER"){
	data.push_back(variable("RUNAFTER", int_type));
      }
      else{
	data.push_back(parse_format(w));
	std::cerr << "\tjoblist variable " << data.back().name << std::endl;
      }
      
    }
    std::cerr << "\treading " << input.size() - 3 << " jobs" << std::endl;
    
    int_data.resize(input.size() - 3);
    real_data.resize(input.size() - 3);
    word_data.resize(input.size() - 3);
    
    for(unsigned int i=2; i<input.size()-1; ++i){
      std::istringstream is(input[i]);
      
      for(unsigned int d = 0; d<data.size(); ++d){
	parse_data(is, data[d], int_data[i-2], real_data[i-2], word_data[i-2]);
      }
    }

    return 0;
  }

  int size() const
  {
    return int_data.size();
  }
  
  std::vector<variable> data;
  std::vector<std::map<std::string, int> > int_data;
  std::vector<std::map<std::string, double> > real_data;
  std::vector<std::map<std::string, std::string> > word_data;
  
};

#endif
