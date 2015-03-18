// utils_DipTraj.h

#ifndef INCLUDED_UTILS_DIPTRAJ
#define INCLUDED_UTILS_DIPTRAJ

#include <string>
#include "../gromos/Exception.h"
#include "../gmath/Vec.h"
#include <vector>

namespace utils{
  class Time;
}

namespace utils{

  class DipTraj_i;
  /**
   * Class DipTraj
   * Defines an instream that can read the GROMOSXX box dipole moment trajectory file.
   *
   * The instream can (currently) handle the TIMESTEP and DIPOLE blocks
   *
   * @class DipTraj
   * @author S. Riniker
   * @sa 
   * @todo Add more block
   */

  /**
   * Class DipData
   * A class to store box dipole moment read in from a restraint trajectory (*.trs)
   *
   * @class DipData
   * @author S. Riniker
   * @ingroup utils
   */
  class DipData{
    public:
      // struct for storing dipole moment and volume
      struct Dip {
        gmath::Vec dipole;
        double vol;
        Dip() : vol(0.0) {}
        Dip(const Dip & dipm) : dipole(dipm.dipole), vol(dipm.vol) {
        }
        Dip & operator=(const Dip & dipm) {
          vol = dipm.vol; dipole = dipm.dipole;
          return *this;
        }
      };
    /**
     * Constructor
     */
    DipData() {}
    /**
     *  DipData copy constructor
     */
    DipData(const DipData &dipdata) {
      m_data = dipdata.m_data;
    }
    /**
     *  DipData deconstructor
     */
    ~DipData() {}
    /**
     * const accessor to dipole moment data
     */
    const std::vector<Dip> & data() const {
      return m_data;
    }
    /**
     * accessor to dipole moment data
     */
    std::vector<Dip> & data() {
      return m_data;
    }
  private:
    std::vector<Dip> m_data;
  };
  
  class DipTraj{
    /**
     * pIMPL
     */
    DipTraj_i *d_this;
    /**
     * skip
     */
    int d_skip;
    /**
     * stride
     */
    int d_stride;
    /**
     * end of file during striding (or skipping)
     */
    int d_stride_eof;

    /**
     * copy constructor
     * not implemented
     */
    DipTraj(const DipTraj&);
    /**
     * operator =
     * not implemented
     */
    DipTraj &operator=(const DipTraj&);
    
  public:
    /**
     * Constructor
     * @param skip   : skip frames
     * @param stride : take only every n-th frame
     */
    DipTraj(int skip=0, int stride=1);
    /**
     * Constructor
     * @param name   : trajectory file
     * @param skip   : skip frames
     * @param stride : take only every n-th frame
     */
    DipTraj(const std::string &name, int skip=0, int stride=1);
    /**
     * Destructor
     */
    ~DipTraj();

    // Methods
    /**
     * open a trajectory file
     * skip and stride continue
     */
    void open(const std::string &name);

    /**
     * read the trajectory
     */
    void read();
    /**
     * close
     */
    void close();
    /**
     * read a jvalue restraints frame
     */
    DipTraj &operator>>(utils::DipData &dipdata);
    /**
     * get the time information from a frame
     */
    DipTraj &operator>>(utils::Time &time);

    // Accessors
    /**
     * title (of the trajectory)
     */
    std::string title()const;
    /**
     * name?
     */
    std::string name()const;
    /**
     * end of file
     */
    bool eof()const;

    /**
     * take every n-th frame
     * starts with frame 0, then n, 2n, 3n, ...
     */
    int stride()const;
    /**
     * skip n frames
     * after skipping, skip will be 0.
     * to skip again, just set it again...
     */
    int skip()const;

    /**
     * set stride
     */
    void stride(int stride);
    /**
     * set skip
     */
    void skip(int skip);

    /**
     * encountered end of file
     * during striding
     */
    bool stride_eof()const;

    //Exceptions
    /**
     * Exception
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what_arg) : gromos::Exception("DipTraj", what_arg){}
    };
  };
}
#endif
