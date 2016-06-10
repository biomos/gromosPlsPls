// pb_Ewald_edir.h


#ifndef INCLUDED_PB_Ewald_edir
#define INCLUDED_PB_Ewald_edir
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif

namespace pb{



class Ewald_edir{
    
bound::Boundary *pbc;
utils::AtomSpecifier atoms;
pb::PB_Parameters ppp;
gcore::Box thebox;

std::vector<std::vector<std::vector< complex<double>  > > >  eir;
//std::vector<std::vector<std::vector< complex<double>  > > >  eir_slow;
bool excluded_interx;

int kmax;
double realcut;
double tolerance;
double box[3];


int kx;
int ky;
int kz;


 public:
  // constructor

  Ewald_edir(bound::Boundary & pbc,utils::AtomSpecifier atoms,
             double realcut, double tolerance,  int KX_in, int KY_in, int KZ_in);// bool exinterx);


   // deconstructor
  ~Ewald_edir(){}



  //methods


     complex<double> scalecomplexnum(       double r,
				            complex<double> c);

     complex<double> multiplycomplexnums(
                                 complex<double> a,
                                 complex<double> b);

     complex<double> multiplycomplexconjugate(
				complex<double> a,
				complex<double> b);


     complex<double> multiplycomplexnums_firstconj(
				complex<double> a,
				complex<double> b);

     complex<double> multiplycomplexnums_bothconj(
				complex<double> a,
				complex<double> b);
	
        void calcenergy(double (& energy)[15]);


        double rspaceEwald(double ewaldcoeff);
	
	
	double calc_ewaldcoeff();
	
	void tabulate_eir(double (& lll)[3]);
	
	double kspaceEwald(
			double ewaldcoeff);


	double kspaceEwald_slow(
			double ewaldcoeff);
        
	double A2timesStildeSquare(
			double ewaldcoeff);

	double minusA3timesStildeSquare(
			double ewaldcoeff);


	double A1timesStildeSquare(
			double ewaldcoeff);


        
        double self_other();
	 double coulomb_non_excluded_atoms();
    double  XIEWcorr();
    //double dip2();
	double minusA1timesSSquare(
			double ewaldcoeff);
	
	void calc_lll(double (& lll)[3]);


        // actually we do not use this one
	double erfcapp(double X);




}; // class
} // namespace


#endif