//proarse program: program to calculate the ARea 'Surface-accessible of' Estimate 

#include <cassert>


#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/physics.h"
#include "../src/fit/PositionUtils.h"
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>



using namespace std;
using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace fit;
using namespace utils;

void heapsort(double* values, int n, int* key);

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "zslice", "sasaatoms", "probe", "traj"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@sasaatoms <atomspecifier: atoms to consider for sasa>\n";
  usage += "\t[@zslice <distance between the Z-slices through the molecule [nm], default: 0.005>]\n";
  usage += "\t[@probe  <probe radius [nm], default: 0.14>]\n";
  usage += "\t@traj <trajectory files>\n";
  
 
    try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    //   get simulation time
  double time=0, dt=1;
  {
    Arguments::const_iterator iter=args.lower_bound("time");
    if(iter!=args.upper_bound("time")){
      time=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("time"))
        dt=atof(iter->second.c_str());
  }

  //try for zslice and probe
    double zslice = 0.005;
    if(args.count("zslice")>0) zslice=atof(args["zslice"].c_str());
    double probe = 0.14;
    if(args.count("probe")>0) probe=atof(args["probe"].c_str());
  
  

  
  //  read topology
  InTopology it(args["topo"]);
  System sys(it.system());
    

  //get sasa atoms
   AtomSpecifier sasaatoms(sys);
    {
       Arguments::const_iterator iter = args.lower_bound("sasaatoms");
       Arguments::const_iterator to = args.upper_bound("sasaatoms");

       for(;iter!=to;iter++){
        string spec=iter->second.c_str();
        sasaatoms.addSpecifier(spec);
       }
    }

    // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);
  // parse gather method
  Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

  //some constants
  double PI = M_PI;
  double twoPI = 2 * PI;        

  //get radii and other things
  AtomSpecifier heavyatoms(sys);
  vector<double> radheavy;
  vector<double> radheavysq;
  int count = 0;
  double rmax = 0;
  int natoms = 0;
  //get all heavy atoms...
  for (int i=0; i < sys.numMolecules(); ++i) {
   natoms += sys.mol(i).numAtoms();
   for (int j=0; j < sys.mol(i).numAtoms(); ++j) {
    if (sys.mol(i).topology().atom(j).mass() != 1.00800) {
     ++count;	    
     heavyatoms.addAtom(i, j);
     double rad = sys.mol(i).topology().atom(j).radius();
     rad +=probe;
     rmax = ((rad > rmax) ? (rad) : (rmax));
     radheavy.push_back(rad);
     radheavysq.push_back(rad*rad);
    }
   }
  }
    // define input coordinate
  InG96 ic;

    
    // loop over all trajectories
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
    ic.open((iter->second).c_str());
      
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      (*pbc.*gathmethod)();
   

        double totSASA = 0;
        int ict = 2000;

	//determine max and min positions using the PositionsUtils
	//using a one-dimensional array for this is not so nice...
        //well, just remember that, MAX position is stored from 0 to 2 (starting at index 0)
        //                          MIN position is stored from 3 to 5



	//this routine is fucked...
        Vec min = fit::PositionUtils::getmincoordinates(&sys, true);
        Vec max = fit::PositionUtils::getmaxcoordinates(&sys, true);

      

        //set up cubicals containing the atoms
	//this is analogous to the hanging spanish moss algorithm from the 70's saturday night live

        double idim = rint((max[0] - min[0])/rmax+0.1);
        if (idim < 3) idim=3;
        double jidim = rint((max[1] - min[1])/rmax+0.1);
        if (jidim < 3) jidim=3;
        jidim = idim * jidim;
        double kjidim = rint((max[2] - min[2])/rmax+0.1);
        if (kjidim < 3) kjidim=3;
        kjidim = jidim * kjidim;


        int kdim = 2000;

        int itab[kdim];
        int natm[natoms][kdim];
        int cube[natoms];
    
        
        

        for (int l=0; l < (int) heavyatoms.size(); ++l) {
         
	  Vec tmp = *heavyatoms.coord(l);
 	  int i = (int) ((tmp[0] - min[0])/rmax + 1.0);
	  int j = (int) ((tmp[1] - min[1])/rmax);
	  int k = (int) ((tmp[2] - min[2])/rmax);

       
         
            // this is the cube index
            int kji = ((int) k) * ((int) jidim) + ((int) j) * ((int) idim) + ((int) i);

	    int n = itab[kji] + 1;

            

            itab[kji] = n;
            natm[n][kji] = l;
            cube[l] = kji;
	    
	}

        
        //go through each atom

	double nzp = rint(0.1/(zslice) + 0.05);

        double dx[2000];
        double dy[2000];
        double dsq[2000]; 
        double d[2000];
        int inov[2000];

        
            
   for (int ir=0; ir < count; ++ir) {

	    int kji = cube[ir];
            int io = 0;
            double area = 0.0;
            Vec tmp =  *heavyatoms.coord(ir);
            double xr = tmp[0];
            double yr = tmp[1];
            double zr = tmp[2];
          
            double rr = radheavy[ir];
            double rrx2 = rr * 2;
            double rrsq = radheavysq[ir];
 
	    double zres = 0;
            double zgrid = 0;

	    int nm = 0;    
	    

	    //loop over cubes
	    // cubeloop:  
               for (int k=-1; k <= 1; ++k) {
                for (int j=-1; j <= 1; ++j) {
                 for (int i=-1; i <= 1; ++i) {
		  int mkji = ((int) kji) + k * ((int) jidim) + j * ((int) idim) + i;
                 
		  if (mkji >= 1) {
		    if(mkji > kjidim) goto breakcube;//break cubeloop; //we break the whole loop...
		    nm=itab[mkji];
		     if(nm >= 1) {
			      //     -- record the atoms in inov that neighbor atom ir
  		      for (int m = 1; m <= nm; ++m) { 
		       int in = natm[m][mkji];
		       if (in != ir) {
		       io += 1;
		       if (io > ict) {
                        ostringstream os;
                        os << "STOP'SOLVA_ERROR: intrsctns > max";
                        throw(gromos::Exception("proarse", os.str()));				     
		       }
		       Vec tmp = *heavyatoms.coord(in);
			dx[io] = xr - tmp[0];
                        dy[io] = yr - tmp[1];
                        dsq[io]= dx[io] * dx[io] + dy[io] * dy[io];
                        d[io] = sqrt(dsq[io]);
                        inov[io] = in;
		       }
		      }
		     }
		  }
		 }
		}
	       }
	    
   breakcube:
   

    bool gotoend = false;

    double b = 0;

    //
    double tf = 0;
    double ti = 0;
      double beta = 0;
    double alpha = 0;
      double trig_test = 0;
      bool sort = false;



    //
 
    double accs[2000];
    double parea = 0;
    double arcsum = 0;
    double t = 0;
    double tt = 0;
    double arcf[2000];
    int karc = 0;
    double rsec2n = 0;
    double rsecn = 0;
        
    if (io >= 1) {
	// z resolution determined
        zres= rrx2/nzp;
        Vec tmp = *heavyatoms.coord(ir);
        zgrid = tmp[2] - rr - zres/2;
    }
    else {
	area = twoPI * rrx2;	
        gotoend = true;
    }

      
    double rsec2r = 0;
    double rsecr = 0;
    double arci[2000];

    
    //     section atom spheres perpendicular to the z axis

    // mainloop: {

      if (gotoend) goto end;//break mainloop; 


     for (int i=0; i < nzp; ++i) {
       // innermainloop: {         

	zgrid = zgrid + zres;

	//     find the radius of the circle of intersection of 
	//     the ir sphere on the current z-plane

	rsec2r = rrsq - (zgrid - zr) * (zgrid - zr);
        rsecr=sqrt(rsec2r);

        
        for (int k=0; k <= karc; ++k) arci[k] = 0.0;
        
         karc = -1;
	 for (int j=1; j <= io; ++j) {
	   //     innerloop: {
	   	               
	    int in = inov[j];
     
	    //find radius of circle locus
            Vec tmp = *heavyatoms.coord(in);   
	    rsec2n = radheavysq[in] - ((zgrid - tmp[2]) * (zgrid - tmp[2]));
           
	    if (rsec2n <= 0.0) goto binnerloop; //break innerloop; 

	      rsecn=sqrt(rsec2n);

              //find intersections of n.circles with ir circles in section

	      if (d[j] >= (rsecr+rsecn))  goto binnerloop; //break innerloop; 
	       	      
	      //do the circles intersect, or is one circle completely inside the other?
	          b = rsecr-rsecn;

		  if (d[j] < abs(b)) {           		      
		    if (b <= 0.0) goto binnermain;//break innermainloop;      	   
	            goto binnerloop; //break innerloop;
		  }

	      //if the circles intersect, find the points of intersection

	        karc += 1;
	      
		 if(karc >= ict) {
                  ostringstream os;
                  os << "STOP'SOLVA_ERROR: karc >= ict";
                  throw(gromos::Exception("proarse", os.str()));		    
		 }

		 
	       //     Initial and final arc endpoints are found for the ir circle intersected
	       //     by a neighboring circle contained in the same plane. The initial endpoint
	       //     of the enclosed arc is stored in arci, and the final arc in arcf
	       //     law of cosines
     
	       trig_test=(dsq[j] + rsec2r-rsec2n)/(2 * d[j] * rsecr);

               if(trig_test >= 1.0)  trig_test=0.99999;
               if(trig_test <= -1.0) trig_test=-0.99999;
               alpha = acos(trig_test);

     
	       //     alpha is the angle between a line containing a point of intersection and
	       //     the reference circle center and the line containing both circle centers
     
               beta = atan2(dy[j],dx[j]) + PI;

	       //     beta is the angle between the line containing both circle centers and the x-axis
     
               ti = beta - alpha;
               tf = beta + alpha;
               if(ti < 0.0) ti += twoPI;
               if(tf > twoPI) tf -= twoPI;
               arci[karc] = ti;
               
               if(tf < ti) { 
        	//if the arc crosses zero, then it is broken into two segments.
	        //the first ends at twoPI and the second begins at zero
               
	       arcf[karc]= twoPI;
               karc += 1;
	       }

               arcf[karc]=tf;	    
	     
	       //	  } // end inner loop
	 binnerloop: continue;
	 } //end for (int j=1; j <= io; ++j) {
     
	 
        //find the accessible surface area for the sphere ir on this section
  
         sort = false;
   
	 //  sumcontrib: {
		  if (karc != -1) { //sort the shit
		    sort = true;
		   }
	           else {	  
                    arcsum= twoPI;		    
                    goto par; //break sumcontrib;
		   }
	      
		//The arc endpoints are sorted on the value of the initial arc endpoint

	        int tag[2000];

        	if (sort) heapsort(arci, karc+1, tag);
		  		                   
		//calculate the accessible area
                 		
		    arcsum=arci[0];
	            t = arcf[tag[0]];		          
                		    
                    if(karc != 0) {
	             for(int k=1; k <= karc; ++k) {                   
		      if(t < arci[k]) arcsum += arci[k]-t;
		      tt=arcf[tag[k]];
		      if(tt > t) t=tt;                   
		     }
		    }		    
		    
                 arcsum += twoPI-t;
                  
     //The area/radius is equal to the accessible arc length x the section thickness.
		 //	 }	     

     par:
              parea=arcsum*zres;
                  
     //Add the accessible area for this atom in this section to the area for this
     //atom for all the section encountered thus far
     
            area += parea;
       	   
	    //   } //end inner main loop
     binnermain: continue;
     } //end for (int i=0; i < nzp; ++i)
     //  } //end mainloop
     
   end:
    
             //scale area to vdw shell

             b=area*rr;
    
             

             accs[ir] = b;
          if (sasaatoms.findAtom(heavyatoms.mol(ir),heavyatoms.atom(ir)) > 0) totSASA += b;
      
   }

   cout << time << ' ' << totSASA << endl;
   time += dt;
    }
    ic.close();
  }
    }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


    void heapsort(double* values, int n, int* key) {

      //this implements a heapsort, which also returns the keys...
    
       int left, right, i, j, hkeys;
       double helper;
       bool ok;
	

  //     initialize index into the original ordering

	      for (i=0; i < n; ++i) key[i] = i;
	      //     perform the heapsort of the input list
	      //		left = values.length/2;
	      //		right = values.length-1;

	      left = n/2;
              right = n-1;

		do {
		 if (left == 0) {
		   helper = values[0];
                   hkeys  = key[0];
		   values[0] = values[right];
                   key[0] = key[right];
		   values[right] = helper;
                   key[right] = hkeys;                                

		   i = 0;
		   right--;
		  } else {
		   left--;
		   i = left;
		  }

		   helper = values[i];
                   hkeys = key[i];
		   ok = false;
		   j = 2 * i;

		   do {
		    if (j < right) {
		     if (values[j] < values[j+1]) j++;
		    }
		     if (helper >= values[j]) {
			ok=true;
		      } else {
		      helper = values[i];
                      hkeys = key[i];
		      values[i] = values[j];
                      key[i] = key[j];
		      values[j] = helper;
                      key[j] = hkeys;
		      i = j;
		      j = 2*i;
		     };
		   } while (!ok && right>=j);

		} while (right >= 1);
       	

    } //end void heapsort
