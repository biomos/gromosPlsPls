/**
 * @file groTime.h
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
   * A class to handle trajectory times. It takes care of the \c \@time command
   * line argument. There are the following options:
   * - no \c \@time argument:
   *  - For programs writing a time series: read the time from the trajectory
   *  - For programs not writing a time series by default: do not write a time
   *    series.
   * - \c \@time argument:
   *  - For programs writing a time series: read the time from the trajectory.
   *  - For programs not writing a time series by default: write a time series
   *    and read the time from the trajectory.
   * - \c \@time \c time \c dt argument: Same behaviour as with the \c \@time
   *   argument but the time is not read from the trajectory but calculated when
   *   a frame is read.
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
     * Constructor for a dummy time object
     */
    Time(): d_t0(0.0), d_dt(1.0), d_current_time(-1.0), d_steps(0), d_read(true), d_do_timeseries(false) {}
    
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
     * get the current steps
     */
    double steps() const {
      return d_steps;
    }
    /**
     * accessor to steps
     */
    double & steps() {
      return d_steps;
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
     * accessor to time series
     */
    bool & doSeries() {
      return d_do_timeseries;
    }
    /**
     * accessor to time series
     */
    bool doSeries() const {
      return d_do_timeseries;
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
     * the current steps
     */
    double d_steps;
    /**
     * read trajectory
     */
    bool d_read;
    /**
     * do a time series: \@time is given but maybe without parameters
     */
    bool d_do_timeseries;
  };
}
#endif
