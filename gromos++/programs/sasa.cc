/**
 * @file sasa.cc
 * Calculates solvent-accessible surface areas for selected atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor sasa
 * @section sasa Calculates solvent-accessible surface areas for selected atoms
 * @author @ref mk
 * @date 21-6-07
 *
 * Program sasa calculates and prints the solvent-accessible surface
 * area (sasa) of all heavy atoms in the solute part of the molecular system.
 * It also calculates the contribution made by a specified set of heavy atoms.
 * The program uses the algorithm of Lee and Richards [J. Mol. Biol., 55, 379-400 (1971)].
 * A spherical probe of given radius is rolled over the surface of the molecule
 * (the size of the probe is typically 0.14~nm for water). The path traced out
 * by its centre gives the accessible surface.
 *
 * In GROMOS, the radii of the heavy atoms are obtained by calculating
 * the minimum energy distance of the interaction between the heavy
 * atom and the first solvent atom. This value is reduced by the specified
 * probe radius to account for the radius of the solvent atom.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;]</td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to consider for sasa&gt; </td></tr>
 * <tr><td> [\@zslice</td><td>&lt;distance between the Z-slices through the molecule (default: 0.005~nm)&gt;] </td></tr>
 * <tr><td> [\@probe</td><td>&lt;probe IAC and radius (default: 4  0.14~nm)&gt;] </td></tr>
 * <tr><td> [\@verbose</td><td>(print summaries)] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory file(s)&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  sasa
    @topo     ex.top
    @pbc      r
    @time     0 1
    @atoms    1:CB
    @zslice   0.005
    @probe    4 0.14
    @verbose
    @traj     ex.tr
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/groTime.h"
#include "../src/utils/AtomicRadii.h"

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

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "time" << "zslice" << "atoms" << "probe" << "traj"
         << "verbose";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <molecular topology file>\n";
  usage += "\t@pbc      <boundary type>\n";
  usage += "\t[@time     <time and dt>]\n";
  usage += "\t@atoms    <atoms to consider for sasa>\n";
  usage += "\t[@zslice  <distance between the Z-slices (default: 0.005)>]\n";
  usage += "\t[@probe   <probe IAC and radius (default: 4 0.14)>]\n";
  usage += "\t[@verbose (print summaries)\n";
  usage += "\t@traj     <trajectory files>\n";
  
  try{
    Arguments args(argc, argv, knowns, usage);

    //   get simulation time
    Time time(args);
    int numFrames=0;
    
    //try for zslice and probe
    double zslice = args.getValue<double>("zslice", false, 0.005);
    vector<double> probearg = args.getValues<double>("probe", 2, false,
          Arguments::Default<double>() << 4 << 0.14);
    int probe_iac = int(probearg[0]);
    double probe = probearg[1];
    
    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    utils::compute_atomic_radii_vdw(probe_iac, probe, sys, it.forceField());
    
    //get sasa atoms
    AtomSpecifier sasaatoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for(;iter!=to;iter++){
        string spec=iter->second.c_str();
        sasaatoms.addSpecifier(spec);
      }
    }

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

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
	if (!sys.mol(i).topology().atom(j).isH()) {
	  ++count;	    
	  heavyatoms.addAtom(i, j);
	  double rad = sys.mol(i).topology().atom(j).radius();
	  rad += probe;
	  rmax = ((rad > rmax) ? (rad) : (rmax));
	  radheavy.push_back(rad);
	  radheavysq.push_back(rad*rad);
	}
      }
    }
    // define input coordinate
    InG96 ic;

    //declare variables, i.e. reserve memory
    int kdim = 2000;
    int ict = 2000;

    vector<int> itab(kdim,0);
    vector<int> inov(kdim,0);
 
    vector<int> empty(kdim, 0);
  
    vector<double> dx(kdim,0.0);
    vector<double> dy(kdim,0.0);
    vector<double> dsq(kdim,0.0);
    vector<double> d(kdim,0.0);
    vector<double> accs(natoms,0.0);
    vector<double> arcf(kdim,0.0);
    vector<double> arci(kdim,0.0);

    vector<int> cube(natoms,0);
    vector<vector<int> > natm;               
    for (int i=0; i < natoms; ++i) natm.push_back(empty);  
	
    // print title
    cout << "#     "
	 << setw(15) << "selected"
	 << setw(10) << "heavy" << endl
	 << "# time"
	 << setw(15) << "atoms"
	 << setw(10) << "atoms" << endl;
    
    

    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){

      // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys >> time;
	(*pbc.*gathmethod)();
   

        double totSASA = 0;
	double totSASA_all = 0;
	
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

        //int itab[kdim];
        //zero out some stuff...

        for (int i=0; i < kdim; ++i) {
	  itab[i] = 0;
	  inov[i] = 0;
	  dx[i] = 0.0;
	  dy[i] = 0.0;
	  dsq[i] = 0.0;
	  d[i] = 0.0;
	}
        
        for (int i=0; i < natoms; ++i) { 
	  cube[i] = 0;
	  natm[i] = empty;  
	}
	
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
            
	for (int ir=0; ir < count; ++ir) {

	  int kji = cube[ir];
	  int io = 0;
	  double area = 0.0;
	  Vec tmp =  heavyatoms.pos(ir);
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
			  throw(gromos::Exception("sasa", os.str()));				     
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

	  std::vector<double> arci(kdim);

	  //
 
	  for (int i=0; i < kdim; ++i) {
	    arcf[i] = 0.0;
	    arci[i] = 0.0;
	  }

    
	  double parea = 0;
	  double arcsum = 0;
	  double t = 0;
	  double tt = 0;
    
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

	    {
	      std::vector<int> tag(kdim);
	      for (int i=0; i < kdim; ++i) tag[i] = 0;

	      if (sort) heapsort(&arci[0], karc+1, &tag[0]);
		  		                   
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
	    }
	    
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
    
	  // add it for averaging
	  accs[ir] += b;
          if (sasaatoms.findAtom(heavyatoms.mol(ir),heavyatoms.atom(ir)) >= 0) totSASA += b;
	  totSASA_all+=b;
	}

	cout.precision(8);
	cout << time;
	cout.precision(5);
	cout << setw(15) << totSASA 
	     << setw(10) << totSASA_all << endl;

	numFrames++;
      }
      ic.close();
    }
    // calculate and print averages
    double totSASA=0.0;
    double totSASA_all=0.0;

    // loop over all heavy atoms
    for(int i=0; i< count; ++i){
      accs[i]/=numFrames;
      totSASA_all+=accs[i];
      if(sasaatoms.findAtom(heavyatoms.mol(i), heavyatoms.atom(i))>=0) totSASA += accs[i];
    }
    cout.precision(5);
    cout << "#\n# ave."
	 << setw(15) << totSASA
	 << setw(10) << totSASA_all << endl;
    if(args.count("verbose") >=0){
	
      cout << "#\n# average contribution per selected heavy atom\n";
      cout << "#\n# "
	   << setw(6) << "atom"
	   << setw(10) << "residue"
	   << setw(5)  << "name" << ' ' 
	   << setw(10) << "SASA" << ' '
	   << setw(10) << "\% selected" << ' '
	   << setw(10) << "\% heavy" << endl;
      
      double sumaccs=0.0, sumper=0.0, sumpera=0.0;
      
      for(int i=0; i< sasaatoms.size(); ++i){
	int index=heavyatoms.findAtom(sasaatoms.mol(i), sasaatoms.atom(i));
	if(index>=0){
	  cout << "# "
	       << setw(6) << sasaatoms.toString(i)
	       << setw(5) << sasaatoms.resnum(i)+1
	       << setw(5) << sasaatoms.resname(i)
	       << setw(5) << sasaatoms.name(i) << ' '
	       << setw(10) << accs[index] << ' '
	       << setw(10) << accs[index]/totSASA * 100.0 << ' '
	       << setw(10) << accs[index]/totSASA_all * 100.0 << ' '
	       << endl;
	  sumaccs+=accs[index];
	  sumper+=accs[index]/totSASA * 100.0;
	  sumpera+=accs[index]/totSASA_all * 100.0;
	}
      }
      cout << "#\n# total                 "
	   << setw(10) << sumaccs << ' '
	   << setw(10) << sumper << ' '
	   << setw(10) << sumpera << ' ' << endl;
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
  //adapted to arrays that start indexing from 0...
  //     initialize index into the original ordering

  for (int i=0; i < n; ++i) key[i] = i;
  //     perform the heapsort of the input list
  //		left = values.length/2;
  //		right = values.length-1;

  int k = n/2 +1;
  int index = n;
  double lists;
  int keys;
    	
  do { 
    		
    if (k > 1) {
      k = k - 1;
      lists = values[k-1];
      keys = key[k-1];
    }
    else {
      lists = values[index-1];
      keys = key[index-1];
      values[index-1] = values[0];
      key[index-1] = key[0];
      index = index - 1;
    }
    if (index <= 1) {
      values[0] = lists;
      key[0] = keys;
      return;
    }
    		
    int i = k;
    int j = k + k;
    do {
    		  
      if (j < index) {
	if (values[j-1] < values[j]) ++j;
      }		
      if (lists < values[j-1]) {
	values[i-1] = values[j-1];
	key[i-1] = key[j-1];
	i = j;
	j = j + j;
      }
      else {
	j = index + 1;
      }
    } while (j <= index);
    			
    values[i-1] = lists;
    key[i-1] = keys;
    			
  } while (n > 1);
       	

} //end void heapsort
