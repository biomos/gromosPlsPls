// pb_FFTGridType.h

#ifndef INCLUDED_PB_FFTGridType
#define INCLUDED_PB_FFTGridType
#ifndef INCLUDED_PB_PB_Parameters
#include "PB_Parameters.h"
#endif



namespace pb{



class FFTGridType{

public:

	/*
	  grid dimensions X 
	 */
	int ngrdx;
	/* 
	  grid dimensions Y 
	 */
	int ngrdy;
	/* 
	  grid dimensions Z 
	 */
	int ngrdz;
	/* 
	  grid dimensions X * grid dimensions Y 
	 */
	int ngrdxy;	
	/*
	  grid dimensions X * grid dimensions Y * grid dimensions Z
	 */
	int ngr3;
	/*
	  X edge of the periodic unit cell
	 */
	double xlen;
	/* 
	  Y edge of the periodic unit cell
	 */
	 double ylen;
	/*
	  Z edge of the periodic unit cell
	 */
	double zlen;
	/*
	  unit cell volume
	 */
	double vol;


        /* the grid center: the middle of the grid*/
        double centerx;
        double centery;
        double centerz;



	/* 
	  X spacing between r-space grid points
	 */               
	 double drx;
	/*
	  Y spacing between r-space grid points
	 */
	double dry;
	/*
	  Z spacing between r-space grid points
	 */
	double drz;
	/*
	  X spacing between k-space grid points
	 */
	double dkx;
	/*
	  Y spacing between k-space grid points
	 */
	double dky;
	/* 
	  Z spacing between k-space grid points
	 */
	double dkz;
	/*
	  number of cubes
	 */
	int ncubes;
	
	
	 //recycles Ewald vacuum field
	
	//boll recycleVACfield = false;
	
	  //writes Ewald vacuum fields to file
	
	//bool writeVACfield = false;
	
	 // reads Ewald vacuum fields to file
	 
	// boolreadVACfield = false;
	// string vacfield1 = "";
	// string vacfield2 = "";
	
	//bool havevacfield = false;



         pb::PB_Parameters ppp;



    //constructor

	
	FFTGridType(int ngrdx, int ngrdy, int ngrdz, double xlen, double ylen, double zlen, int ncubes);
        FFTGridType(){}
    // deconstructor
       ~FFTGridType(){}

        //methods

	void dumpparameters();

        
}; // class
} // namespace


#endif
