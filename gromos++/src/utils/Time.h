/**
 * @file Time.h
 * io of TIMESTEP and \@time
 */

#ifndef INCLUDED_UTILS_TIME
#define INCLUDED_UTILS_TIME

#include <iostream>

namespace args {
  class Arguments;
}

namespace gio {
  class InG96;
}

namespace utils{
  /**
   * Class Time
   * A class to handle trajectory times.
   *
   * @class Time
   * @author A. Eichenberger, N. Schmid
   * @ingroup utils
   */
  class Time{
  public:
    /**
     * Constructor
     * @param args the command line arguments
     */
    Time(const args::Arguments & args);
    
    /**
     * get the current time
     */
    double time() const {
      return d_current_time;
    }
    /**
     * accessor to time
     */
    double & time() {
      return d_current_time;
    }
    
    /**
     * accessor to timestep
     */
    double dt() const {
      return d_dt;
    }
    /**
     * accessor to timestep
     */
    double & dt() {
      return d_dt;
    }
    
    /**
     * accessor to starting time
     */
    double start_time() const {
      return d_t0;
    }
    /**
     * accessor to starting time
     */
    double & start_time() {
      return d_t0;
    }
    /**
     * accessor to read
     */
    bool & read() {
      return d_read;
    }
    /**
     * accessor to read
     */
    bool read() const {
      return d_read;
    }
    
    /**
     * print the time to a stream
     * @param out the stream
     */
    void print(std::ostream & out) const;
    
    /**
     * Prints the time to a stream
     */
    friend std::ostream & operator<<(std::ostream &os, Time const & t);
  protected:
    /**
     * the time argument t0
     */
    double d_t0;
    /**
     * the time argument dt
     */
    double d_dt;
    /**
     * the current time
     */
    double d_current_time;
    /**
     * read trajectory
     */
    bool d_read;
  };
}
#endif
