/**
 * @file Time.cc
 * implements utils::Time
 */

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "../args/Arguments.h"
#include "Time.h"

namespace utils {
  Time::Time(const args::Arguments & args) :
  d_t0(0.0), d_dt(1.0), d_current_time(0.0), d_read(true) {
    // get the time command line argument
    args::Arguments::const_iterator iter = args.lower_bound("time");
    if (iter != args.upper_bound("time")) {
      // argument supplied: do not read
      read() = false;
      std::istringstream in(iter->second);
      if (!(in >> start_time()))
        throw gromos::Exception("time", "@time: starting time is not numeric.");
      ++iter;
    }
    if (iter != args.upper_bound("time")) {
      std::istringstream in(iter->second);
      if (!(in >> dt()))
        throw gromos::Exception("time", "@time: timestep is not numeric.");
    }

    if (read() == false) {
      // the first frame has time start_time(). so when we call the operator
      // to read the time it will be increased to start_time() + dt(). That's
      // why we have to start at start_time() - dt()
      time() = start_time() - dt();
    }
  }

  void Time::print(std::ostream & out) const {
    // write the time
    out.setf(std::ios::fixed, std::ios::floatfield);
    out << std::setw(15) << std::setprecision(9) << time();
  }
  
  std::ostream & operator<<(std::ostream &os, Time const & t) {
    t.print(os);
    return os;
  }
}




