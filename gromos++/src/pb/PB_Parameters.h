// pb_PB_Parameters.h

#ifndef INCLUDED_PB_PB_Parameters
#define INCLUDED_PB_PB_Parameters

namespace pb{



class PB_Parameters{

public:
 //grid center
 //double gridcenterX ;
 //double gridcenterY;
// double gridcenterZ ;
 
 //grid points
 //int GPX;
 //int GPY;
 //int GPZ;
 
 //grid spacing
 //double gridspacing;
 
 //permittivity
 double epssolute; 
// double epssolvent;

 double FPEPSI;

 double PI;

 double tiny_real;
 
 double tinfoil_boundary;
 double adjusted_boundary;
 double vacuum_boundary;

 double epssolvent;

 // convergence limit for change in deltaG, in FFT Iterator
 double ddg_conv;

 //hydrogen radius [nm]
 double hradius;

 // max number of surface points (FFT)
 int nmax_surfpoints;

 // default charge shaping function (FFT)
 int default_chshape; 

 // widths double alpha1, double alpha2, and alias vector numbers, int nalias1, int nalias2 for the FFT algorithm

 double alpha1;
 double alpha2;
 int nalias1;
 int nalias2;

 // relaxation parameter for FFT iterative algorithm
 double FFTlambda;

 // default thread number for FFT
 int threadnum;

 double kappa;


 int debugvar; // 0 or 1: no or yes debugging

 double convergence_fd;
 double convergence_fft;


 double tole;
 double rcut;
 int kmax;

 double quadr;


 // xi ewald
 double xiew;

public:
  // constructor
  PB_Parameters(double epssolvent);
  PB_Parameters();

   // deconstructor
  ~PB_Parameters(){}

// methods

 double getPI();
 double getFPEPSI();
 double getEpssolute();
 double getTiny_real();
 double get_rf_boundary_eps();
 double get_sc_boundary_eps();
 double get_ls_boundary_eps();
 double get_hradius();
 int get_nmax_surfpoints();
 int get_default_chshape();

 int get_nalias1();
 int get_nalias2();
 double get_alpha1();
 double get_alpha2();
 double get_FFTlambda();

 int get_threadnum();
 double getKappa();

 int get_debugvar();

 double get_convergence_fd();
 double get_convergence_fft();


  double get_rcut();
  double get_tole();
  int get_kmax();

  double get_quadr();
  double get_xiew();

  
  // private:
  void initParams(double epssolvent);



}; //class
}//namespace
#endif