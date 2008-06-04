// mk_script.h

enum filetype{unknownfile, inputfile, topofile, coordfile, refposfile, 
	      posresspecfile, disresfile, pttopofile, dihresfile, jvaluefile,
	      ledihfile, outputfile, outtrxfile, outtrvfile, outtrefile, outtrgfile, 
	      scriptfile, outbaefile, outbagfile};
int numFiletypes=19;
typedef map<string, filetype>::value_type FT;
const FT filetypes[] ={FT("", unknownfile),
		       FT("input", inputfile),
		       FT("topo", topofile),
		       FT("coord", coordfile),
		       FT("refpos", refposfile),
		       FT("posresspec", posresspecfile),
		       FT("disres", disresfile),
		       FT("pttopo", pttopofile),
		       FT("dihres", dihresfile),
		       FT("jvalue", jvaluefile),
		       FT("ledih", ledihfile),
		       FT("output", outputfile),
		       FT("outtrx", outtrxfile),
		       FT("outtrv", outtrvfile),
		       FT("outtre", outtrefile),
		       FT("outtrg", outtrgfile),
		       FT("outbae", outbaefile),
		       FT("outbag", outbagfile),
		       FT("script", scriptfile)
};
static map<string,filetype> FILETYPE(filetypes,filetypes+numFiletypes);

enum blocktype {unknown, systemblock, startblock, initialiseblock,
                minimiseblock, energyminblock, stochasticblock, stochdynblock,
                readtrajblock, consistencycheckblock, stepblock, boundaryblock,
                boundcondblock, multicellblock, submoleculesblock, tcoupleblock,
                thermostatblock, multibathblock, pcoupleblock, barostatblock,
                virialblock, pressurescaleblock, centreofmassblock,
                overalltransrotblock, comtransrotblock,
		        printblock, printoutblock, writeblock, writetrajblock,
		        ewarnblock, debugblock,
                shakeblock, geomconstraintblock, constraintblock, forceblock,
                covalentformblock, plistblock, neighbourlistblock, pairlistblock,
                plist03block, nonbondedblock, longrangeblock, cgrainblock,
                posrestblock, positionresblock, distrestblock, distanceresblock,
                diherestblock, dihedralresblock, jvalblock, jvalueresblock,
                localelevationblock, localelevblock, rottransblock, perturbblock,
                perturbationblock, lambdasblock, perturb03block, umbrellablock,
                perscaleblock, replicablock, innerloopblock, fourdimblock,
                gromos96compatblock, pathintblock, integrateblock, randomnumbersblock};

typedef map<string, blocktype>::value_type BT;
int numBlocktypes = 67;
const BT blocktypes[] ={BT("",unknown),
			BT("SYSTEM",systemblock),
			BT("START",startblock),
            BT("INITIALISE",initialiseblock),
			BT("MINIMISE",minimiseblock),
            BT("ENERGYMIN",energyminblock),
            BT("STOCHASTIC",stochasticblock),
            BT("STOCHDYN",stochdynblock),
            BT("READTRAJ",readtrajblock),
            BT("CONSISTENCYCHECK",consistencycheckblock),
            BT("STEP",stepblock),
			BT("BOUNDARY",boundaryblock),
            BT("BOUNDCOND",boundcondblock),
            BT("MULTICELL",multicellblock),
			BT("SUBMOLECULES",submoleculesblock),
			BT("TCOUPLE",tcoupleblock),
            BT("THERMOSTAT",thermostatblock),
            BT("MULTIBATH",multibathblock),
			BT("PCOUPLE",pcoupleblock),
            BT("BAROSTAT",barostatblock),
            BT("VIRIAL",virialblock),
            BT("PRESSURESCALE",pressurescaleblock),
			BT("CENTREOFMASS",centreofmassblock),
            BT("OVERALLTRANSROT",overalltransrotblock),
            BT("COMTRANSROT",comtransrotblock),
			BT("PRINT",printblock),
            BT("PRINTOUT",printoutblock),
			BT("WRITE",writeblock),
            BT("WRITETRAJ",writetrajblock),
            BT("EWARN",ewarnblock),
            BT("DEBUG",debugblock),
            BT("SHAKE",shakeblock),
            BT("GEOMCONSTRAINTS",geomconstraintblock),
            BT("CONSTRAINT",constraintblock),
			BT("FORCE",forceblock),
            BT("COVALENTFORM",covalentformblock),
			BT("PLIST",plistblock),
            BT("NEIGHBOURLIST",neighbourlistblock),
            BT("PAIRLIST",pairlistblock),
			BT("PLIST03",plist03block),
            BT("NONBONDED",nonbondedblock),
			BT("LONGRANGE",longrangeblock),
            BT("CGRAIN",cgrainblock),
			BT("POSREST",posrestblock),
            BT("POSITIONRES",positionresblock),
            BT("DISTREST",distrestblock),
            BT("DISTANCERES",distanceresblock),
            BT("DIHEREST",diherestblock),
            BT("DIHEDRALRES",dihedralresblock),
            BT("J-VAL",jvalblock),
            BT("JVALUERES",jvalueresblock),
            BT("LOCALELEVATION",localelevationblock),
            BT("LOCALELEV",localelevblock),
            BT("ROTTRANS",rottransblock),
			BT("PERTURB",perturbblock),
            BT("PERTURBATION",perturbationblock),
            BT("LAMBDAS", lambdasblock),
			BT("PERTURB03",perturb03block),
            BT("UMBRELLA",umbrellablock),
            BT("PERSCALE",perscaleblock),
            BT("REPLICA",replicablock),
            BT("INNERLOOP",innerloopblock),
            BT("FOURDIM",fourdimblock),
            BT("GROMOS96COMPAT",gromos96compatblock),
            BT("PATHINT",pathintblock),
            BT("INTEGRATE",integrateblock),
			BT("RANDOMNUMBERS",randomnumbersblock),
};
static map<string,blocktype> BLOCKTYPE(blocktypes,blocktypes+numBlocktypes);

enum templateelement{unknowntemplate, systemtemplate, numbertemplate,
		     oldnumbertemplate,  
		     start_timetemplate, end_timetemplate, queuetemplate};
typedef map<string, templateelement>::value_type TE;
const TE templateelements[]={TE("", unknowntemplate),
			     TE("system", systemtemplate),
			     TE("number", numbertemplate),
			     TE("oldnumber", oldnumbertemplate),
			     TE("start_time", start_timetemplate),
			     TE("end_time", end_timetemplate),
			     TE("queue", queuetemplate)
};
static map<string,templateelement> TEMPLATE(templateelements, 
					    templateelements+7);

//BLOCKDEFINITIONS
class isystem{
public:
  int found, npm, nsm;
  isystem(){found=0;}
};

class istart{
public:
  int found, ntx,init,ntx0;
  double ig,tempi,heat,boltz;
  istart(){found=0;}
};

class iinitialise{
public:
  int found, ntivel, ntishk, ntinht, ntinhb;
  int ntishi, ntirtc, nticom, ntisti;
  double ig, tempi;
  iinitialise(){found=0;}
};

class iminimise{
public:
  int found, ntem, ncyc, nmin;
  double dele, dx0, dxm;
  iminimise(){found=0;}
};

class ienergymin{
public:
  int found, ntem, ncyc, nmin;
  double dele, dx0, dxm, flim;
  ienergymin(){found=0;}
};

class istochastic{
public:
  int found, ntfr, nsfr, nbref;
  double rcutf, cfric, tempsd;
  istochastic(){found=0;}
};

class istochdyn{
public:
  int found, ntsd, ntfr, nsfr, nbref;
  double rcutf, cfric, tempsd;
  istochdyn(){found=0;}
};

class ireadtraj{
public:
  int found, ntrd, ntrn, ntrb, ntshk;
  ireadtraj(){found=0;}
};

class iconsistencycheck{
public:
  int found, ntchk, ntckf, ntckv, ntckt;
  int ntcke, ntckr, ntckl;
  double fdckf, fdckv, fdckl;
  vector<int> nckf;
  iconsistencycheck(){found=0;}
};

class istep{
public:
  int found, nstlim;
  double t,dt;
  istep(){found=0;}
};

class iboundary{
public:
  int found, ntb,nrdbox;
  double box[3],beta;
  iboundary(){found=0;}
};

class iboundcond{
public:
  int found, ntb, ndfmin;
  iboundcond(){found=0;}
};

class imulticell{
public:
  int found, ntm, ncella, ncellb, ncellc;
  double tolpx, tolpv, tolpf, tolpfw;
  imulticell(){found=0;}
};

class isubmolecules{
public:
  int found;
  vector<int> nsp;
  isubmolecules(){found=0;}
};

class itcouple{
public:
  int found, ntt[3];
  double temp0[3], taut[3];
  itcouple(){found=0;}
};

class ithermostat{
public:
  int found, ntt, ntvar;
  class tbath{
  public:
    double tembth;
    vector<double> taubth;
  };
  vector<tbath> baths; 
  class dofgroup{
  public:
    int ntcpl, nttyp, ntcns, ntgt;
    vector<int> ntgtg;
  };
  vector<dofgroup> dofgroups;
  ithermostat(){found=0;}
};

class imultibath{
public:
  int found, num, nbaths, dofset;
  string algorithm;
  vector<double> temp0, tau;
  vector<int> last, combath, irbath;
  imultibath(){found=0; num=-1;}
};

class ipcouple{
public:
  int found, ntp;
  double pres0,comp, taup;
  ipcouple(){found=0;}
};

class ibarostat{
public:
  int found, ntp, npvar, npcpl[6];
  double comp;
  class pbath{
  public:
    double prsbth;
    vector<double> taubba; 
  };
  vector<pbath> pbaths;
  ibarostat(){found=0;} 
};

class ivirial{
public:
  int found, ntv, ntvg;
  ivirial(){found=0;}
};

class ipressurescale{
public:
  int found;
  string couple, scale, virial;
  double comp, taup;
  double pres0[3][3];
  ipressurescale(){found=0;}
};

class icentreofmass{
public:
  int found, ndfmin,ntcm,nscm;
  icentreofmass(){found=0;}
};

class ioveralltransrot{
public:
  int found, ncmtr, ncmro;
  double cmamx, cmamy, cmamz;
  ioveralltransrot(){found=0;}
};

class icomtransrot{
public:
  int found, nscm;
  icomtransrot(){found=0;}
};

class iprint{
public:
  int found, ntpr, ntpl, ntpp;
  iprint(){found=0;}
};

class iprintout{
public:
  int found, ntpr, ntpp;
  iprintout(){found=0;}
};

class iwrite{
public:
  int found, ntwx,ntwse,ntwv,ntwe,ntwg,ntba,ntpw;
  iwrite(){found=0; ntwx=0; ntwse=0; ntwv=0; ntwe=0; ntwg=0; ntba=-1; ntpw=0;}
};

class iwritetraj{
public:
  int found, ntwx, ntwse, ntwv, ntwf, ntwe, ntwg, ntwb;
  iwritetraj(){found=0; ntwx=0; ntwse=0; ntwv=0; ntwf=0; ntwe=0; ntwg=0; ntwb=0;}
};

class iewarn{
public:
  int found;
  double maxener;
  iewarn(){found=0;}
};

class idebug{
public:
  int found;
  class routine{
  public:
     string piider;
     int iiideo;    
  };
  vector<routine> routines;
  idebug(){found=0;}
};

class ishake{
public:
  int found, ntc;
  double tol;
  ishake(){found=0;}
};

class igeomconstraint{
public:
  int found, ntcph, ntcpn, ntcs;
  double shktol;
  igeomconstraint(){found=0;}
}; 

class iconstraint{
public:
  int found, ntc;
  double ntcp0[3], ntcs0[3]; 
  string ntcp, ntcs;
  iconstraint(){found=0;}
};

class iforce{
public:
  int found, ntf[10];
  vector<int> nre;
  iforce(){found=0;}
};

class icovalentform{
public:
  int found, ntbbh, ntbah, ntbdn;
  icovalentform(){found=0;}
};

class iplist{
public:
  int found, ntnb, nsnb;
  double rcutp, rcutl;
  iplist(){found=0;}
};

class ineighbourlist{
public:
  int found, nplsr, nmesr, ntrsr, npllr, nmelr, ntrlr, nplgx, nplgy, nplgz;
  double rlstp, rcutp, rlsts, rcuts, rlstl, rcutl;
  ineighbourlist(){found=0;}
};

class ipairlist{
public:
  int found, nsnb;
  double rcutp, rcutl;
  string algorithm, type, size;
  ipairlist(){found=0;}
};

class iplist03{
public:
  int found, nsnb;
  double rcutp, rcutl, grds;
  bool chargegroup, grid, grda;
  iplist03(){found=0;}
};

class inonbonded{
public:
  int found, nlrele, nshape, na2clc, nkx, nky, nkz, ngx, ngy, ngz;
  int nasord, nfdord, nalias, nqeval, nrdgrd, nwrgrd, nlrlj;
  double appak, rcrf, epsrf, ashape, tola2, epsls, kcut, nspord, faccur, slvdns;
  inonbonded(){found=0;}
};

class ilongrange{
public:
  int found;
  double epsrf, appak, rcrf;
  ilongrange(){found=0;}
};

class icgrain{
public:
  int found, ntcgran;
  double eps;
  icgrain(){found=0;}
};

class iposrest{
public:
  int found, ntr, nrdrx;
  double cho;
  iposrest(){found=0;}
};

class ipositionres{
public:
  int found, ntpor, ntporb, ntpors;
  double cpor;
  ipositionres(){found=0;}
};

class idistrest{
public:
  int found, ntdr, nrddr;
  double cdis, dr0, taudr;
  idistrest(){found=0;}
};

class idistanceres{
public:
  int found, ntdir, ntdira;
  double cdir, dir0, taudir;
  idistanceres(){found=0;}
};

class idiherest{
public:
  int found, ntdlr;
  double cdlr;
  idiherest(){found=0;}
};

class idihedralres{
public:
  string ntdlr;
  int found, philin;
  double cdlr;
  idihedralres(){found=0;}
};

class ijval{
public:
  int found, ntjr, ntjrh, nrdjr;
  double cjr, taujr;
  ijval(){found=0;}
};

class ijvalueres{
public:
  int found, ntjvra, le, ngrid;
  string ntjvr;
  double cjvr, taujvr, delta;
  ijvalueres(){found=0;}
};

class ilocalelevation{
public:
  int found, ntle, nrdle;
  double cwle;
  ilocalelevation(){found=0;}
};

class ilocalelev{
public:
  int found, ntles, ntlesa;
  double cles;
  ilocalelev(){found=0;}
};

class irottrans{
public:
  int found, rtc, rtclast;
  irottrans(){found=0;}
};

class iperturb{
public:
  int found, ntg,nrdgl,nlam,mmu;
  double rlam, dlamt, rmu, dmut,alphlj,alphc;
  iperturb(){found=0;}
};

class iperturbation{
public:
  int found, ntg, nrdgl, nlam, nscale;
  double rlam, dlamt, alphlj, alphc;
  iperturbation(){found=0;}
};

class ilambdas{
public:
  int found, ntil;
  class lambint{
  public:
    string ntli;
    int nilg1, nilg2;
    double ali, bli, cli, dli, eli;
  };
  vector<lambint> lambints;
  ilambdas(){found=0;}
};

class iperturb03{
public:
  int found, ntg, nlam, scaling;
  double rlam, dlamt,alphlj,alphc;
  iperturb03(){found=0;}
};

class iumbrella{
public:
  int found, ntus;
  double uscst1, uscst2, usref1, usref2;
  iumbrella(){found=0;}
};

class iperscale{
public:
  int found, t, read;
  double kdih, kj, diff, ratio;
  string restype;
  iperscale(){found=0;}
};

class ireplica{
public:
  int found, lrescale, nretrial, nrequil, nrejob, nrewrt;
  vector<double> ret, relam, rets; 
  ireplica(){found=0; nrewrt=0;}
};

class iinnerloop{
public:
  int found, spec;
  iinnerloop(){found=0;}
};

class ifourdim{
public:
  int found, nt4dim, ndfmi4, nt4xi, nt4xo, ntt4, ntcw4d;
  int ntf4[6];
  double cw4da, temp4i, temp04, taut4, cw4db, temp0b;
  ifourdim(){found=0;}
};

class igromos96compat{
public:
  int found, ntnb96, ntr96, ntp96, ntg96;
  igromos96compat(){found=0;}
};

class ipathint{
public:
  int found, ntpi;
  ipathint(){found=0;}
};

class iintegrate{
public:
  int found, nint;
  iintegrate(){found=0;}
};

class irandomnumbers{
public:
  int found, ntrng, ntgsl;
  irandomnumbers(){found=0;}
};

class iunknown
{
public:
  string name;
  string content;
  iunknown(string n):name(n){}
};


class input{
public:
  isystem system;
  istart start;
  iinitialise initialise;
  iminimise minimise;
  ienergymin energymin;
  istochastic stochastic;
  istochdyn stochdyn;
  ireadtraj readtraj;
  iconsistencycheck consistencycheck;
  istep step;
  iboundary boundary;
  iboundcond boundcond;
  imulticell multicell;
  isubmolecules submolecules;
  itcouple tcouple;
  ithermostat thermostat;
  imultibath multibath;
  ipcouple pcouple;
  ibarostat barostat;
  ivirial virial;
  ipressurescale pressurescale;
  icentreofmass centreofmass;
  ioveralltransrot overalltransrot;
  icomtransrot comtransrot;
  iprint print;
  iprintout printout;
  iwrite write;
  iwritetraj writetraj;
  iewarn ewarn;
  idebug debug;
  ishake shake;
  igeomconstraint geomconstraint;
  iconstraint constraint;
  iforce force;
  icovalentform covalentform;
  iplist plist;
  ineighbourlist neighbourlist;
  ipairlist pairlist;
  iplist03 plist03;
  inonbonded nonbonded;
  ilongrange longrange;
  icgrain cgrain;
  iposrest posrest;
  ipositionres positionres;
  idistrest distrest;
  idistanceres distanceres;
  idiherest diherest;
  idihedralres dihedralres;
  ijval jval;
  ijvalueres jvalueres;
  ilocalelevation localelevation;
  ilocalelev localelev;
  irottrans rottrans;
  iperturb perturb;
  iperturbation perturbation;
  ilambdas lambdas;
  iperturb03 perturb03;
  iumbrella umbrella;
  iperscale perscale;
  ireplica replica;
  iinnerloop innerloop;
  ifourdim fourdim;
  igromos96compat gromos96compat;
  ipathint pathint;
  iintegrate integrate;
  irandomnumbers randomnumbers;
  
  vector<iunknown> unknown;
}; 

class fileInfo{
public:
  double box[3];
  vector<string> blocks;
  vector<int> blockslength;
};

//INSTREAM
istringstream &operator>>(istringstream &is, isystem &s){
  string e;
  s.found=1;
  is >> s.npm >> s.nsm >> e;
  return is;
}
istringstream &operator>>(istringstream &is, istart &s){
  string e;
  s.found=1;
  is >> s.ntx >> s.init >> s.ig >> s.tempi >> s.heat >> s.ntx0 >> s.boltz >> e;
  return is;
}
istringstream &operator>>(istringstream &is, iinitialise &s){
  string e;
  s.found=1;
  is >> s.ntivel >> s.ntishk >> s.ntinht >> s.ntinhb;
  is >> s.ntishi >> s.ntirtc >> s.nticom;
  is >> s.ntisti >> s.ig >> s.tempi >> e;
  return is;
}
istringstream &operator>>(istringstream &is, iminimise &s){
  string e;
  s.found=1;
  is >> s.ntem >> s.ncyc >> s.dele >> s.dx0 >> s.dxm;
  if (!(is >> s.nmin)){
    s.nmin = 0;
  }
  else
    is >> e;
	
  return is;
}
istringstream &operator>>(istringstream &is, ienergymin &s){
  string e;
  s.found=1;
  is >> s.ntem >> s.ncyc >> s.dele >> s.dx0 >> s.dxm >> s.nmin >> s.flim >> e;
  return is;
}
istringstream &operator>>(istringstream &is, istochastic &s){
  string e;
  s.found=1;
  is >> s.ntfr >> s.nsfr >> s.nbref >> s.rcutf >> s.cfric >> s.tempsd >> e;
  return is;
}
istringstream &operator>>(istringstream &is, istochdyn &s){
  string e;
  s.found=1;
  is >> s.ntsd >> s.ntfr >> s.nsfr >> s.nbref >> s.rcutf >> s.cfric >> s.tempsd >> e;
  return is;
}
istringstream &operator>>(istringstream &is, ireadtraj &s){
  string e;
  s.found=1;
  is >> s.ntrd >> s.ntrn >> s.ntrb >> s.ntshk >> e;
  return is;
}
istringstream &operator>>(istringstream &is, iconsistencycheck &s){
  string e;
  s.found=1;
  is >> s.ntchk >> s.ntckf >> s.fdckf >> s.ntckv >> s.fdckv >> s.ntckt
     >> s.ntcke >> s.ntckr >> s.ntckl >> s.fdckl;
  int nackf, nckf;
  is >> nackf;
  for(int i=0; i<nackf; i++) {is >> nckf; s.nckf.push_back(nckf);}
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is, istep &s){
  string e;
  s.found=1;
  is >> s.nstlim >> s.t >> s.dt >> e;
  return is;
}
istringstream &operator>>(istringstream &is, iboundary &s){
  string e;
  s.found=1;
  is >> s.ntb >> s.box[0] >> s.box[1] >> s.box[2] >> s.beta >> s.nrdbox >> e;
  return is;
}
istringstream &operator>>(istringstream &is, iboundcond &s){
  string e;
  s.found=1;
  is >> s.ntb >> s.ndfmin >> e;
  return is;
}
istringstream &operator>>(istringstream &is, imulticell &s){
  string e;
  s.found=1;
  is >> s.ntm >> s.ncella >> s.ncellb >> s.ncellc
     >> s.tolpx >> s.tolpv >> s.tolpf >> s.tolpfw >> e;
  return is;
}
istringstream &operator>>(istringstream &is, isubmolecules &s){
  string e;
  s.found=1;
  int nspm,nsp;
  is >> nspm;
  for(int i=0; i<nspm; i++){is >> nsp; s.nsp.push_back(nsp);}
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is,itcouple &s){
  string e;
  s.found=1;
  for(int i=0;i<3;i++) is >> s.ntt[i] >> s.temp0[i] >> s.taut[i];
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ithermostat &s){
  string e;
  int dum, ntbth, ntgrp;
  double tau;
  s.found=1;
  is >> s.ntt >> s.ntvar >> ntbth >> ntgrp;
  for(int i=0;i<ntbth;++i){
    class ithermostat::tbath b;
    is >> dum >> b.tembth;
    for(int j=0;j<s.ntvar; ++j){
       is >> tau;
       b.taubth.push_back(tau);
    }
    s.baths.push_back(b);
  }
  for(int i=0;i<ntgrp;++i){
    class ithermostat::dofgroup d;
    is >> d.ntcpl >> d.nttyp >> d.ntcns >> d.ntgt;
    for(int j=0;j<d.ntgt; ++j){
      is >> dum;
      d.ntgtg.push_back(dum); 
    }
    s.dofgroups.push_back(d);
  }
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is,imultibath &s){
  string e, alg;
  double temp0, tau;
  int last, combath, irbath;
  s.found=1;
  is >> alg;
  std::transform(alg.begin(), alg.end(), alg.begin(), ::tolower);
  s.algorithm=alg;
  if(alg=="nose-hoover-chains"){
    is >> s.num;
  }

  is >> s.nbaths;
  for(int i=0;i<s.nbaths;++i){
    is >> temp0;
    s.temp0.push_back(temp0);
  }
  for(int i=0;i<s.nbaths;++i){
    is >> tau;
    s.tau.push_back(tau);
  }

  is >> s.dofset;
  for(int i=0;i<s.dofset;++i){
    is >> last;
    s.last.push_back(last);
  }
  for(int i=0;i<s.dofset;++i){
    is >> combath;
    s.combath.push_back(combath);
  }
  for(int i=0;i<s.dofset;++i){
    is >> irbath;
    s.irbath.push_back(irbath);
  }
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ipcouple &s){
  string e;
  s.found=1;
  is >> s.ntp >> s.pres0 >> s.comp >> s.taup >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ibarostat &s){
  string e;
  int dum, npbth;
  double taup;
  s.found=1;
  is >> s.ntp >> s.npvar >> npbth >> s.comp;
  for(int i=0;i<npbth;++i){
    class ibarostat::pbath p;
    is >> dum >> p.prsbth;
    for(int j=0;j<s.npvar; ++j){
       is >> taup;
       p.taubba.push_back(taup);
    }
    s.pbaths.push_back(p);
  }
  is >> s.npcpl[0] >> s.npcpl[1] >> s.npcpl[2] >> s.npcpl[3] >> s.npcpl[4] >> s.npcpl[5] >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ivirial &s){
  string e;
  s.found=1;
  is >> s.ntv >> s.ntvg >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ipressurescale &s){
  string e, coup, sca, vir;
  s.found=1;

  is >> coup;
  std::transform(coup.begin(), coup.end(), coup.begin(), ::tolower);
  s.couple=coup;

  is >> sca;
  std::transform(sca.begin(), sca.end(), sca.begin(), ::tolower);
  s.scale=sca;

  is >>  s.comp >> s.taup >> vir;
  std::transform(vir.begin(), vir.end(), vir.begin(), ::tolower);
  s.virial=vir;
  
  is >> s.pres0[0][0] >> s.pres0[0][1] >> s.pres0[0][2]
     >> s.pres0[1][0] >> s.pres0[1][1] >> s.pres0[1][2]
     >> s.pres0[2][0] >> s.pres0[2][1] >> s.pres0[2][2] >> e;

  return is;
}
istringstream &operator>>(istringstream &is,icentreofmass &s){
  string e;
  s.found=1;
  is >> s.ndfmin >> s.ntcm >> s.nscm >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ioveralltransrot &s){
  string e;
  s.found=1;
  is >> s.ncmtr >> s.ncmro >> s.cmamx >> s.cmamy >> s.cmamz >> e;
  return is;
}
istringstream &operator>>(istringstream &is,icomtransrot &s){
  string e;
  s.found=1;
  is >> s.nscm >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iprint &s){
  string e;
  s.found=1;
  is >> s.ntpr >> s.ntpl >> s.ntpp >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iprintout &s){
  string e;
  s.found=1;
  is >> s.ntpr >> s.ntpp >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iwrite &s){
  string e;
  s.found=1;
  is >> s.ntwx >> s.ntwse >> s.ntwv >> s.ntwe >> s.ntwg >> s.ntba;

  if (!(is >> s.ntpw)){
    // std::cout << "could not read ntpw, assume g96" << std::endl;
    s.ntpw = s.ntba;
    s.ntba = -1;
  }
  else{
    // std::cout << "could read ntba, ntpw! => gXX" << std::endl;
    is >> e;
  }
    
  return is;
}
istringstream &operator>>(istringstream &is,iwritetraj &s){
  string e;
  s.found=1;
  is >> s.ntwx >> s.ntwse >> s.ntwv >> s.ntwf >> s.ntwe >> s.ntwg >> s.ntwb >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iewarn &s){
  string e;
  s.found=1;
  is >> s.maxener >> e;
  return is;
}
istringstream &operator>>(istringstream &is,idebug &s){
  string e, piid;
  int nrd;
  s.found=1;
  is >> nrd;
  for(int i=0; i<nrd; ++i){
    class idebug::routine r;
    is >> piid;
    // std::transform(piid.begin(), piid.end(), piid.begin(), ::tolower);
    r.piider=piid;
    is >> r.iiideo;
    s.routines.push_back(r);
  }
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ishake &s){
  string e;
  s.found=1;
  is >> s.ntc >> s.tol >> e;
  return is;
}
istringstream &operator>>(istringstream &is,igeomconstraint &s){
  string e;
  s.found=1;
  is >> s.ntcph >> s.ntcpn >> s.ntcs >> s.shktol >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iconstraint &s){
  string e, cp, cs;
  s.found=1;
  is >> s.ntc >> cp;
  std::transform(cp.begin(), cp.end(), cp.begin(), ::tolower);
  s.ntcp=cp;
  is >> s.ntcp0[0];
  if(s.ntcp=="flexshake"){
    is >> s.ntcp0[1] >> s.ntcp0[2];
  }
  is >> cs;
  std::transform(cs.begin(), cs.end(), cs.begin(), ::tolower);
  s.ntcs=cs;
  is >> s.ntcs0[0];
  if(s.ntcs=="flexshake"){
    is >> s.ntcs0[1] >> s.ntcs0[2];
  }
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iforce &s){
  string e;
  s.found=1;
  int negr, nre;
  for(int i=0; i<10; i++) is >> s.ntf[i];
  is >> negr;
  for(int i=0; i<negr; i++) {is >> nre; s.nre.push_back(nre);}
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is,icovalentform &s){
  string e;
  s.found=1;
  is >> s.ntbbh >> s.ntbah >> s.ntbdn >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iplist &s){
  string e;
  s.found=1;
  is >> s.ntnb >> s.nsnb >> s.rcutp >> s.rcutl >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ineighbourlist &s){
  string e;
  s.found=1;
  is >> s.nplsr >> s.nmesr >> s.ntrsr >> s.rlstp >> s.rcutp;
  is >> s.npllr >> s.nmelr >> s.ntrlr >> s.rlsts >> s.rcuts >> s.rlstl >> s.rcutl;
  is >> s.nplgx >> s.nplgy >> s.nplgz >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ipairlist &s){
  string e, alg, siz, typ;
  s.found=1;
  is >> alg;
  std::transform(alg.begin(), alg.end(), alg.begin(), ::tolower);
  s.algorithm = alg;
  is >> s.nsnb >> s.rcutp >> s.rcutl >> siz;
  std::transform(siz.begin(), siz.end(), siz.begin(), ::tolower);
  s.size = siz;
  is >> typ;
  std::transform(typ.begin(), typ.end(), typ.begin(), ::tolower);
  s.type = typ;
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is, iplist03 &s){
  string e, alg, grdsz, ctoff;
  s.found = 1;
  is >> alg >> s.nsnb >> s.rcutp >> s.rcutl >> grdsz >> ctoff >> e;

  std::transform(alg.begin(), alg.end(), alg.begin(), ::tolower);
  std::transform(grdsz.begin(), grdsz.end(), grdsz.begin(), ::tolower);
  std::transform(ctoff.begin(), ctoff.end(), ctoff.begin(), ::tolower);

  if (ctoff == "chargegroup") s.chargegroup = true;
  else if (ctoff == "atomic") s.chargegroup = false;
  if (grdsz == "auto") s.grda = true;
  else{
    std::istringstream is(grdsz);
    is >> s.grds;
  }
  if (alg == "standard") s.grid = false;
  else if (alg == "grid") s.grid = true;
  else throw gromos::Exception("input", "wrong pairlist algorithm");
  return is;
}
istringstream &operator>>(istringstream &is,inonbonded &s){
  string e;
  s.found=1;
  is >> s.nlrele;
  is >> s.appak >> s.rcrf >> s.epsrf;
  is >> s.nshape >> s.ashape >> s.na2clc >> s.tola2 >> s.epsls;
  is >> s.nkx >> s.nky >> s.nkz >> s.kcut;
  is >> s.ngx >> s.ngy >> s.ngz >> s.nasord >> s.nfdord >> s.nalias >> s.nspord;
  is >> s.nqeval >> s.faccur >> s.nrdgrd >> s.nwrgrd;
  is >> s.nlrlj >> s.slvdns >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ilongrange &s){
  string e;
  s.found=1;
  is >> s.epsrf >> s.appak >> s.rcrf >> e;
  return is;
}
istringstream &operator>>(istringstream &is,icgrain &s){
  string e;
  s.found=1;
  is >> s.ntcgran >> s.eps >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iposrest &s){
    string e;
    s.found=1;
    is >> s.ntr >> s.cho >> s.nrdrx >> e;
    return is;
}
istringstream &operator>>(istringstream &is,ipositionres &s){
  string e;
  s.found=1;
  is >> s.ntpor >> s.ntporb >> s.ntpors >> s.cpor >> e;
  return is;
}
istringstream &operator>>(istringstream &is,idistrest &s){
  string e;
  s.found=1;
  is >> s.ntdr >> s.cdis >> s.dr0 >> s.taudr >> s.nrddr >> e;
  return is;
}
istringstream &operator>>(istringstream &is,idistanceres &s){
  string e;
  s.found=1;
  is >> s.ntdir >> s.ntdira >> s.cdir >> s.dir0 >> s.taudir >> e;
  return is;
}
istringstream &operator>>(istringstream &is,idiherest &s){
  string e;
  s.found=1;
  is >> s.ntdlr >> s.cdlr >> e;
  return is;
}
istringstream &operator>>(istringstream &is,idihedralres &s){
  string e, dlr;
  s.found=1;
  is >> dlr;
  std::transform(dlr.begin(), dlr.end(), dlr.begin(), ::tolower);
  s.ntdlr = dlr;
  is >> s.cdlr >> s.philin >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ijval &s){
  string e;
  s.found=1;
  is >> s.ntjr >> s.ntjrh >> s.cjr >> s.taujr >> s.nrdjr >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ijvalueres &s){
  string e, jvr;
  s.found=1;
  is >> jvr;
  std::transform(jvr.begin(), jvr.end(), jvr.begin(), ::tolower);
  s.ntjvr = jvr;
  is >> s.ntjvra >> s.cjvr >> s.taujvr >> s.le >> s.ngrid >> s.delta >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ilocalelevation &s){
  string e;
  s.found=1;
  is >> s.ntle >> s.cwle >> s.nrdle >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ilocalelev &s){
  string e;
  s.found=1;
  is >> s.ntles >> s.ntlesa >> s.cles >> e;
  return is;
}
istringstream &operator>>(istringstream &is,irottrans &s){
  string e;
  s.found=1;
  is >> s.rtc >> s.rtclast >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iperturb &s){
  string e;
  s.found=1;
  is >> s.ntg >> s.nrdgl >> s.rlam >> s.dlamt >> s.rmu >> s.dmut
     >> s.alphlj >> s.alphc >> s.nlam >> s.mmu >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iperturbation &s){
  string e;
  s.found=1;
  is >> s.ntg >> s.nrdgl >> s.rlam >> s.dlamt
     >> s.alphlj >> s.alphc >> s.nlam >> s.nscale >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ilambdas &s){
  string e, dum;
  int ntil;
  s.found=1;
  is >> ntil;
  for(int i=0; i<ntil; ++i){
    is >> dum;
    std::transform(dum.begin(), dum.end(), dum.begin(), ::tolower);
    class ilambdas::lambint l;
    l.ntli=dum; 
    is >> l.nilg1 >> l.nilg2 >> l.ali >> l.bli >> l.cli >> l.dli >> l.eli;
    s.lambints.push_back(l);
  }
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iperturb03 &s){
  string e;
  s.found=1;
  is >> s.ntg >> s.rlam >> s.dlamt >> s.nlam
     >> s.alphlj >> s.alphc
     >> s.scaling >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iumbrella &s){
  string e;
  s.found=1;
  is >> s.ntus >> s.uscst1 >> s.uscst2 >> s.usref1 >> s.usref2 >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iperscale &s){
  string e, dum;
  s.found=1;
  is >> dum;
  std::transform(dum.begin(), dum.end(), dum.begin(), ::tolower); 
  s.restype=dum;
  is >> s.kdih >> s.kj >> s.t >> s.diff >> s.ratio >> s.read >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ireplica &s){
  string e;
  int nret, nrelam;
  double dum;
  s.found=1;
  is >> nret;
  for(int i=0; i<nret; ++i){
    is >> dum;
    s.ret.push_back(dum);
  }
  is >> nrelam;
  for(int i=0; i<nrelam; ++i){
    is >> dum;
    s.relam.push_back(dum);
    is >> dum;
    s.rets.push_back(dum);
  }
  is >> s.nretrial >> s.nrequil >> s.nrejob >> s.nrewrt >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iinnerloop &s){
  string e;
  s.found=1;
  is >> s.spec >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ifourdim &s){
  string e;
  s.found=1;
  is >> s.nt4dim >> s.cw4da >> s.temp4i >> s.ndfmi4 >> s.nt4xi >> s.nt4xo;
  is >> s.ntt4 >> s.temp04 >> s.taut4;
  is >> s.ntcw4d >> s.cw4db >> s.temp0b;
  is >> s.ntf4[0] >> s.ntf4[1] >> s.ntf4[2] >> s.ntf4[3] >> s.ntf4[4] >> s.ntf4[5] >> e;
  return is;
}
istringstream &operator>>(istringstream &is,igromos96compat &s){
  string e;
  s.found=1;
  is >> s.ntnb96 >> s.ntr96 >> s.ntp96 >> s.ntg96 >> e;
  return is;
}
istringstream &operator>>(istringstream &is,ipathint &s){
  string e;
  s.found=1;
  is >> s.ntpi >> e;
  return is;
}
istringstream &operator>>(istringstream &is,iintegrate &s){
  string e;
  s.found=1;
  is >> s.nint >> e;
  return is;
}

istringstream &operator>>(istringstream &is, irandomnumbers &s){
  string e;
  s.found=1;
  is >> s.ntrng >> s.ntgsl >> e;
  return is;
}

istringstream &operator>>(istringstream &is, iunknown &s){
  string e=is.str();
  s.content=e.substr(0,e.find("END"));
  return is;
}


Ginstream &operator>>(Ginstream &is,input &gin){
  vector<string> buffer;
  while (!is.stream().eof()){
    is.getblock(buffer);
    if(buffer.size()){
      
    
      string bufferstring;
      gio::concatenate(buffer.begin()+1, buffer.end()-1, bufferstring);
      istringstream bfstream(bufferstring);
      switch(BLOCKTYPE[buffer[0]]){
	case systemblock:           bfstream >> gin.system;           break;
	case startblock:            bfstream >> gin.start;            break;
    case initialiseblock:       bfstream >> gin.initialise;       break;
	case minimiseblock:         bfstream >> gin.minimise;         break;
    case energyminblock:        bfstream >> gin.energymin;        break;
    case stochasticblock:       bfstream >> gin.stochastic;       break;
    case stochdynblock:         bfstream >> gin.stochdyn;         break;
    case readtrajblock:         bfstream >> gin.readtraj;         break;
    case consistencycheckblock: bfstream >> gin.consistencycheck; break;
	case stepblock:             bfstream >> gin.step;             break;
	case boundaryblock:         bfstream >> gin.boundary;         break;
    case boundcondblock:        bfstream >> gin.boundcond;        break;
    case multicellblock:        bfstream >> gin.multicell;        break;
	case submoleculesblock:     bfstream >> gin.submolecules;     break;
	case tcoupleblock:          bfstream >> gin.tcouple;          break;
    case thermostatblock:       bfstream >> gin.thermostat;       break;
    case multibathblock:        bfstream >> gin.multibath;        break;
	case pcoupleblock:          bfstream >> gin.pcouple;          break;
    case barostatblock:         bfstream >> gin.barostat;         break;
    case virialblock:           bfstream >> gin.virial;           break;
    case pressurescaleblock:    bfstream >> gin.pressurescale;    break;
	case centreofmassblock:     bfstream >> gin.centreofmass;     break;
    case overalltransrotblock:  bfstream >> gin.overalltransrot;  break;
    case comtransrotblock:      bfstream >> gin.comtransrot;      break;
	case printblock:            bfstream >> gin.print;            break;
    case printoutblock:         bfstream >> gin.printout;         break;
	case writeblock:            bfstream >> gin.write;            break;
    case writetrajblock:        bfstream >> gin.writetraj;        break;
    case ewarnblock:            bfstream >> gin.ewarn;            break;
    case debugblock:            bfstream >> gin.debug;            break;
	case shakeblock:            bfstream >> gin.shake;            break;
    case geomconstraintblock:   bfstream >> gin.geomconstraint;   break;
    case constraintblock:       bfstream >> gin.constraint;       break;
	case forceblock:            bfstream >> gin.force;            break;
    case covalentformblock:     bfstream >> gin.covalentform;     break;
	case plistblock:            bfstream >> gin.plist;            break;
    case neighbourlistblock:    bfstream >> gin.neighbourlist;    break;
    case pairlistblock:         bfstream >> gin.pairlist;         break;
	case plist03block:          bfstream >> gin.plist03;          break;
    case nonbondedblock:        bfstream >> gin.nonbonded;        break;
	case longrangeblock:        bfstream >> gin.longrange;        break;
    case cgrainblock:           bfstream >> gin.cgrain;           break;
	case posrestblock:          bfstream >> gin.posrest;          break;
    case positionresblock:      bfstream >> gin.positionres;      break;
	case distrestblock:         bfstream >> gin.distrest;         break;
    case distanceresblock:      bfstream >> gin.distanceres;      break;
    case diherestblock:         bfstream >> gin.diherest;         break;
    case dihedralresblock:      bfstream >> gin.dihedralres;      break;
    case jvalblock:             bfstream >> gin.jval;             break;
    case jvalueresblock:        bfstream >> gin.jvalueres;        break;
    case localelevationblock:   bfstream >> gin.localelevation;   break;
    case localelevblock:        bfstream >> gin.localelev;        break;
    case rottransblock:         bfstream >> gin.rottrans;         break;
	case perturbblock:          bfstream >> gin.perturb;          break;
    case perturbationblock:     bfstream >> gin.perturbation;     break;
    case lambdasblock:          bfstream >> gin.lambdas;          break;
	case perturb03block:        bfstream >> gin.perturb03;        break;
    case umbrellablock:         bfstream >> gin.umbrella;         break;
    case perscaleblock:         bfstream >> gin.perscale;         break;
    case replicablock:          bfstream >> gin.replica;          break;
    case innerloopblock:        bfstream >> gin.innerloop;        break;
    case fourdimblock:          bfstream >> gin.fourdim;          break;
    case gromos96compatblock:   bfstream >> gin.gromos96compat;   break;
    case pathintblock:          bfstream >> gin.pathint;          break;
    case integrateblock:        bfstream >> gin.integrate;        break;
	case randomnumbersblock:    bfstream >> gin.randomnumbers;    break;
	  
	case unknown:
	  iunknown newblock(buffer[0]);
	  bfstream >> newblock;
	  gin.unknown.push_back(newblock);
	  cout << "Don't know anything about block " << buffer[0]
	       << ". Just storing data.\n";
      }
    }
    
  }
  return is;
}

Ginstream &operator>>(Ginstream &is,fileInfo &s){

  string e;
  string first;
  vector<string> buffer;
  is.getline(first);
    
  while(!is.stream().eof()){
    is.getblock(buffer);
    s.blocks.push_back(first);
    if (first=="BOX") {
      istringstream iss(buffer[0]);
      iss >> s.box[0] >> s.box[1] >> s.box[2] >> e;
    }
    if (first=="TRICLINICBOX"){
      double dummy1, dummy2;
      istringstream iss(buffer[1]);
      iss >> s.box[0];
      iss.str(buffer[2]);
      iss >> dummy1 >> s.box[1];
      iss.str(buffer[3]);
      iss >> dummy1 >> dummy2 >> s.box[2];
    }
      
    s.blockslength.push_back(buffer.size()-1);
    is.getline(first);
  }
  return is;
}

// TEMPLATE handling of (output) filenames

class filename
{
  vector<string> d_parts;
  double d_time, d_dt;
  int d_start;
  string d_system;
  string d_queue;
  string d_template;
  
public:
  filename();
  filename(string s, double t, double dt, int start=1, string q=""){
    d_system=s;
    d_time=t;
    d_dt=dt;
    d_start=start;
    d_queue=q;
  };
  void setInfo(string s, double t, double dt, int start=1, string q=""){
    d_system=s;
    d_time=t;
    d_dt=dt;
    d_start=start;
    d_queue=q;
  };
  
  void setTemplate(string s);
  string temp(){ return d_template;};
  string name(int number);
};

void filename::setTemplate(string s)
{
  d_template=s;
  d_parts.clear();
  string::size_type iter;
  
  string sub;
  iter=s.find('%');
  while(iter !=string::npos){
    sub=s.substr(0,iter);
    s=s.substr(iter+1,s.size()-iter-1);
    iter=s.find('%');
    d_parts.push_back(sub);
    
  }
  d_parts.push_back(s);
}

string filename::name(int number)
{
  ostringstream os;
  for(unsigned int i=0; i<d_parts.size();i++){
    if(i%2){
      switch(TEMPLATE[d_parts[i]]){
	case systemtemplate:     os << d_system;               break;
	case numbertemplate:     os << d_start + number;       break;
	case oldnumbertemplate:  os << d_start + number -1 ;   break;
	case start_timetemplate: os << d_time+number*d_dt; break;
	case end_timetemplate:   os << d_time+(number+1)*d_dt;     break;
        case queuetemplate:      os << d_queue;                break;
	case unknowntemplate:
	  cout << "Do not know how to handle " << d_parts[i] 
	       << " in template. Just printing the words." << endl;
	  os << d_parts[i];
	  break;
      }
    }
    else os << d_parts[i];
  }
  return os.str();
}

// Jobinfo
class jobinfo
{
public:
  map<string, string> param;
  string dir;
  int prev_id;
};


// Writing out of an input file
ostream &operator<<(ostream &os, input &gin)
{
  // SYSTEM (g96, promd, md++)
  if(gin.system.found)
    os << "SYSTEM\n"
       << "#      NPM      NSM\n"
       << setw(10) << gin.system.npm
       << setw(9)  << gin.system.nsm
       << "\nEND\n";
  // MINIMISE (g96): only write if minimise is on
  if(gin.minimise.found && gin.minimise.ntem){
    os << "MINIMISE\n"
       << "#    NTEM    NCYC    DELE    DX0    DXM\n"
       << setw(10) << gin.minimise.ntem
       << setw(10) << gin.minimise.ncyc
       << setw(10) << gin.minimise.dele
       << setw(10) << gin.minimise.dx0
       << setw(10) << gin.minimise.dxm;
    if (gin.minimise.nmin)
      os << setw(10) << gin.minimise.nmin;
    os << "\nEND\n";
  }
  // ENERGYMIN (promd, md++): only write if NTEM != 0
  if(gin.energymin.found && gin.energymin.ntem){
    os << "ENERGYMIN\n"
       << "#    NTEM    NCYC    DELE    DX0    DXM   NMIN   FLIM\n"
       << setw(10) << gin.energymin.ntem
       << setw(10) << gin.energymin.ncyc
       << setw(10) << gin.energymin.dele
       << setw(10) << gin.energymin.dx0
       << setw(10) << gin.energymin.dxm
       << setw(10) << gin.energymin.nmin
       << setw(10) << gin.energymin.flim
       << "\nEND\n";
  }
  // STOCHASTIC (g96)
  if(gin.stochastic.found){
    os << "STOCHASTIC\n"
       << "#     NTFR      NSFR     NBREF     RCUTF     CFRIC    TEMPSD\n"
       << setw(10) << gin.stochastic.ntfr
       << setw(10) << gin.stochastic.nsfr
       << setw(10) << gin.stochastic.nbref
       << setw(10) << gin.stochastic.rcutf
       << setw(10) << gin.stochastic.cfric
       << setw(10) << gin.stochastic.tempsd
       << "\nEND\n";
  }
  // STOCHDYN (promd, md++)
  if(gin.stochdyn.found){
    os << "STOCHDYN\n"
       << "#     NTSD      NTFR      NSFR     NBREF     RCUTF     CFRIC    TEMPSD\n"
       << setw(10) << gin.stochdyn.ntsd
       << setw(10) << gin.stochdyn.ntfr
       << setw(10) << gin.stochdyn.nsfr
       << setw(10) << gin.stochdyn.nbref
       << setw(10) << gin.stochdyn.rcutf
       << setw(10) << gin.stochdyn.cfric
       << setw(10) << gin.stochdyn.tempsd
       << "\nEND\n";
  }
  // START (g96)
  if(gin.start.found){
    os << "START\n"
       << "# values for INIT -- starting procedures in RUNMD.f\n"
       << "#  INIT\n"
       << "#     shake X\n"
       << "#            shake V\n"
       << "# centre of mass motion removal if NTCM=1\n"
       << "#     1   yes    yes   yes\n"
       << "#     2    no    yes   yes\n"
       << "#     3    no     no   yes\n"
       << "#     4    no     no    no\n"
       << "#      NTX     INIT      IG     TEMPI      HEAT  NTXO   BOLTZ\n"
       << setw(10) << gin.start.ntx
       << setw(9)  << gin.start.init
       << setw(8)  << gin.start.ig
       << setw(10)  << gin.start.tempi
       << setw(10) << gin.start.heat
       << setw(6)  << gin.start.ntx0
       << setw(12)  << gin.start.boltz
       << "\nEND\n";
  }
  // INITIALISE (promd, md++)
  if(gin.initialise.found){
    os << "INITIALISE\n"
       << "# Default values for NTI values: 0\n"
       << "#   NTIVEL    NTISHK    NTINHT    NTINHB\n"
       << setw(10) << gin.initialise.ntivel
       << setw(10) << gin.initialise.ntishk
       << setw(10) << gin.initialise.ntinht
       << setw(10) << gin.initialise.ntinhb
       << "\n"
       << "#   NTISHI    NTIRTC    NTICOM\n"
       << setw(10) << gin.initialise.ntishi
       << setw(10) << gin.initialise.ntirtc
       << setw(10) << gin.initialise.nticom
       << "\n"
       << "#   NTISTI\n"
       << setw(10) << gin.initialise.ntisti
       << "\n"
       << "#       IG     TEMPI\n"
       << setw(10) << gin.initialise.ig
       << setw(10) << gin.initialise.tempi
       << "\nEND\n";
  }
  // READTRAJ (promd, md++)
  if(gin.readtraj.found){
    os << "READTRAJ\n"
       << "#     NTRD      NTRN      NTRB     NTSHK\n"
       << setw(10) << gin.readtraj.ntrd
       << setw(10) << gin.readtraj.ntrn
       << setw(10) << gin.readtraj.ntrb
       << setw(10) << gin.readtraj.ntshk
       << "\nEND\n";
  }
  // CONSISTENCYCHECK (promd)
  if(gin.consistencycheck.found){
    os << "CONSISTENCYCHECK\n"
       << "#    NTCHK     NTCKF     FDCKF     NTCKV     FDCKV\n"
       << setw(10) << gin.consistencycheck.ntchk
       << setw(10) << gin.consistencycheck.ntckf
       << setw(10) << gin.consistencycheck.fdckf
       << setw(10) << gin.consistencycheck.ntckv
       << setw(10) << gin.consistencycheck.fdckv
       << "\n"
       << "#    NTCKT     NTCKE     NTCKR\n"
       << setw(10) << gin.consistencycheck.ntckt
       << setw(10) << gin.consistencycheck.ntcke
       << setw(10) << gin.consistencycheck.ntckr
       << "\n"
       << "#    NTCKL     FDCKL\n"
       << setw(10) << gin.consistencycheck.ntckl
       << setw(10) << gin.consistencycheck.fdckl
       << "\n"
       << "#    NACKF\n"
       << setw(10) << gin.consistencycheck.nckf.size()
       << "\n"
       << "# NCKF(1...NACKF)\n";
    for(unsigned int i=0; i<gin.consistencycheck.nckf.size(); ++i){
       os << setw(10) << gin.consistencycheck.nckf[i];
    }
    os << "\nEND\n";
  }
  // STEP (promd, md++, g96)
  if(gin.step.found){
    os << "STEP\n"
       << "#   NSTLIM         T        DT\n"
       << setw(10) << gin.step.nstlim
       << setw(10) << gin.step.t
       << setw(10) << gin.step.dt
       << "\nEND\n";
  }
  // BOUNDARY (g96)
  if(gin.boundary.found){
    os << "BOUNDARY\n"
       << "#      NTB    BOX(1)    BOX(2)    BOX(3)      BETA  NRDBOX\n"
       << setw(10) << gin.boundary.ntb
       << setw(10) << gin.boundary.box[0]
       << setw(10) << gin.boundary.box[1] 
       << setw(10) << gin.boundary.box[2]
       << setw(10) << gin.boundary.beta
       << setw(8) << gin.boundary.nrdbox
       << "\nEND\n";
  }
  // BOUNDCOND (promd, md++)
  if(gin.boundcond.found){
    os << "BOUNDCOND\n"
       << "#      NTB    NDFMIN\n"
       << setw(10) << gin.boundcond.ntb
       << setw(10) << gin.boundcond.ndfmin
       << "\nEND\n";
  }
  // MULTICELL (promd, md++)
  if(gin.multicell.found){
    os << "MULTICELL\n"
       << "#      NTM    NCELLA    NCELLB    NCELLC\n"
       << setw(10) << gin.multicell.ntm
       << setw(10) << gin.multicell.ncella
       << setw(10) << gin.multicell.ncellb
       << setw(10) << gin.multicell.ncellc
       << "\n"
       << "#     TOLPX    TOLPV     TOLPF    TOLPFW\n"
       << setw(10) << gin.multicell.tolpx
       << setw(10) << gin.multicell.tolpv
       << setw(10) << gin.multicell.tolpf
       << setw(10) << gin.multicell.tolpfw
       << "\nEND\n";
  }
  // SUBMOLECULES (g96)
  if(gin.submolecules.found){
    os << "SUBMOLECULES\n"
       << "#     NSPM  NSP(1.. NSPM)\n"
       << setw(10) << gin.submolecules.nsp.size() << "\n";
    for(unsigned int i=0; i< gin.submolecules.nsp.size(); i++){
      os << setw(6) << gin.submolecules.nsp[i];
      if((i+1)%10==0) os << "\n";
    }
    os << "\nEND\n";
  }
  // TCOUPLE (g96)
  if(gin.tcouple.found){
    os << "TCOUPLE\n"
       << "#      NTT     TEMP0      TAUT\n"
       << setw(10) << gin.tcouple.ntt[0] 
       << setw(10) << gin.tcouple.temp0[0] 
       << setw(10) << gin.tcouple.taut[0] << "\n"
       << setw(10) << gin.tcouple.ntt[1] 
       << setw(10) << gin.tcouple.temp0[1] 
       << setw(10) << gin.tcouple.taut[1] << "\n"
       << setw(10) << gin.tcouple.ntt[2] 
       << setw(10) << gin.tcouple.temp0[2] 
       << setw(10) << gin.tcouple.taut[2]
       << "\nEND\n";
  }
  // THERMOSTAT (promd)
  if(gin.thermostat.found){
    os << "THERMOSTAT\n"
       << "#      NTT     NTVAR     NTBTH     NTGRP\n"
       << setw(10) << gin.thermostat.ntt
       << setw(10) << gin.thermostat.ntvar
       << setw(10) << gin.thermostat.baths.size()
       << setw(10) << gin.thermostat.dofgroups.size()
       << "\n"
       << "# I = 1 ... NTBTH\n"
       << "# I  TEMBTH(I)  TAUBTH(I,1...NTVAR)\n";
    for(unsigned int i=0; i<gin.thermostat.baths.size(); ++i){
      os << setw(3)  << i+1
         << setw(11) << gin.thermostat.baths[i].tembth;
      for(unsigned int j=0; j<gin.thermostat.baths[i].taubth.size(); ++j){
        os << setw(7) << gin.thermostat.baths[i].taubth[j];
      }
      os << "\n";
    }
    os << "# J = 1 ... NTGRP\n"
       << "# NTCPL(J)  NTTYP(J)  NTCNS(J)   NTGT(J)  NTGTG(J,1...NTGT(J))\n";
    for(unsigned int j=0; j<gin.thermostat.dofgroups.size(); ++j){
      os << setw(10) << gin.thermostat.dofgroups[j].ntcpl
         << setw(10) << gin.thermostat.dofgroups[j].nttyp
         << setw(10) << gin.thermostat.dofgroups[j].ntcns
         << setw(10) << gin.thermostat.dofgroups[j].ntgt;
      if(gin.thermostat.dofgroups[j].ntgt>0){
        for(unsigned int k=0; k<gin.thermostat.dofgroups[j].ntgtg.size(); ++k){
          os << setw(7) << gin.thermostat.dofgroups[j].ntgtg[k];
        }
      }
      os << "\n"; 
    }
    os << "END\n";
  }
  // MULTIBATH (md++)
  if(gin.multibath.found){
    os << "MULTIBATH\n"
       << "# ALGORITHM:\n"
       << "#      weak-coupling:      use weak-coupling scheme\n"
       << "#      nose-hoover:        use Nose Hoover scheme\n"
       << "#      nose-hoover-chains: use Nose Hoover chains scheme\n"
       << "# NUM: number of chains in Nose Hoover chains scheme\n"
       << "#      !! only specify NUM when needed !!\n"
       << "# NBATHS: number of temperature baths to couple to\n";
    if(gin.multibath.algorithm=="nose-hoover-chains"){
       os << "#          ALGORITHM     NUM\n"
          << setw(20) << gin.multibath.algorithm
          << setw(8) << gin.multibath.num;
    }
    else{
       os << "#          ALGORITHM\n"
          << setw(20) << gin.multibath.algorithm;
    }
    os << "\n#  NBATHS\n"
       << setw(10) << gin.multibath.nbaths
       << "\n";
    os << "# TEMP0(1 ... NBATHS)\n";
    for(int i=0; i<gin.multibath.nbaths; ++i){
      os << setw(10) << gin.multibath.temp0[i];
    }
    os << "\n";
    os << "# TAU(1 ... NBATHS)\n";
    for(int i=0; i<gin.multibath.nbaths; ++i){
      os << setw(10) << gin.multibath.tau[i];
    }
    os << "\n";
    os << "#   DOFSET: number of distiguishable sets of d.o.f.\n";
    os << setw(10) << gin.multibath.dofset << "\n";
    os << "# LAST(1 ... DOFSET)\n";
    for(int i=0; i<gin.multibath.dofset; ++i){
      os << setw(10) << gin.multibath.last[i];
    }
    os << "\n";
    os << "# COM-BATH(1 ... DOFSET)\n";
    for(int i=0; i<gin.multibath.dofset; ++i){
      os << setw(10) << gin.multibath.combath[i];
    }
    os << "\n";
    os << "# IR-BATH(1 ... DOFSET)\n";
    for(int i=0; i<gin.multibath.dofset; ++i){
      os << setw(10) << gin.multibath.irbath[i];
    }
    os << "\nEND\n";
  }
  // PCOUPLE (g96)
  if(gin.pcouple.found){
    os << "PCOUPLE\n"
       << "#      NTP     PRES0      COMP      TAUP\n"
       << setw(10) << gin.pcouple.ntp
       << setw(10) << gin.pcouple.pres0
       << setw(10) << gin.pcouple.comp
       << setw(10) << gin.pcouple.taup
       << "\nEND\n";
  }
  // BAROSTAT (promd)
  if(gin.barostat.found){
    os << "BAROSTAT\n"
       << "#      NTP      NPVAR     NPBTH      COMP\n"
       << setw(10) << gin.barostat.ntp
       << setw(10) << gin.barostat.npvar
       << setw(10) << gin.barostat.pbaths.size()
       << setw(10) << gin.barostat.comp
       << "\n"
       << "# I = 1 ... NPBTH\n"
       << "# I  PRSBTH(I)  TAUBBA(I,1...NPVAR)\n";
    for(unsigned int i=0; i<gin.barostat.pbaths.size(); ++i){
      os << setw(3)  << i+1
         << setw(11) << gin.barostat.pbaths[i].prsbth;
      for(unsigned int j=0; j<gin.barostat.pbaths[i].taubba.size(); ++j){
        os << setw(7) << gin.barostat.pbaths[i].taubba[j];
      }
      os << "\n";
    }
    os << "#  NPCPL(1...6)\n"
       << setw(6) << gin.barostat.npcpl[0]
       << setw(6) << gin.barostat.npcpl[1]
       << setw(6) << gin.barostat.npcpl[2]
       << setw(6) << gin.barostat.npcpl[3]
       << setw(6) << gin.barostat.npcpl[4]
       << setw(6) << gin.barostat.npcpl[5]
       << "\nEND\n";
  }
  // VIRIAL (promd)
  if(gin.virial.found){
    os << "VIRIAL\n"
       << "#      NTV       NTVG\n"
       << setw(10) << gin.virial.ntv
       << setw(10) << gin.virial.ntvg
       << "\nEND\n";
  }
  // PRESSURESCALE
  if(gin.pressurescale.found){
    os << "PRESSURESCALE\n"
       << "# COUPLE   SCALE    COMP    TAUP  VIRIAL\n"
       << setw(8) << gin.pressurescale.couple
       << setw(8) << gin.pressurescale.scale
       << setw(8) << gin.pressurescale.comp
       << setw(8) << gin.pressurescale.taup
       << setw(8) << gin.pressurescale.virial
       << "\n# PRES0(1...3,1...3)"
       << setw(8) << gin.pressurescale.pres0[0][0]
       << setw(8) << gin.pressurescale.pres0[0][1]
       << setw(8) << gin.pressurescale.pres0[0][2]
       << "\n"
       << setw(8) << gin.pressurescale.pres0[1][0]
       << setw(8) << gin.pressurescale.pres0[1][1]
       << setw(8) << gin.pressurescale.pres0[1][2]
       << "\n"
       << setw(8) << gin.pressurescale.pres0[2][0]
       << setw(8) << gin.pressurescale.pres0[2][1]
       << setw(8) << gin.pressurescale.pres0[2][2]
       << "\nEND\n";
  } 
  // CENTREOFMASS (g96)
  if(gin.centreofmass.found){
    os << "CENTREOFMASS\n"
       << "#   NDFMIN      NTCM      NSCM\n"
       << setw(10) << gin.centreofmass.ndfmin
       << setw(10) << gin.centreofmass.ntcm
       << setw(10) << gin.centreofmass.nscm
       << "\nEND\n";
  }
  // OVERALLTRANSROT (promd)
  if(gin.overalltransrot.found){
    os << "OVERALLTRANSROT\n"
       << "#    NCMTR     NCMRO     CMAMX     CMAMY     CMAMZ\n"
       << setw(10) << gin.overalltransrot.ncmtr
       << setw(10) << gin.overalltransrot.ncmro
       << setw(10) << gin.overalltransrot.cmamx
       << setw(10) << gin.overalltransrot.cmamy
       << setw(10) << gin.overalltransrot.cmamz
       << "\nEND\n";
  };
  // COMTRANSROT (md++)
  if(gin.comtransrot.found){
    os << "COMTRANSROT\n"
       << "#     NSCM\n"
       << setw(10) << gin.comtransrot.nscm
       << "\nEND\n";
  }
  // PRINT (g96)
  if(gin.print.found){
    os << "PRINT\n"
       << "#NTPR: print out energies, etc. every NTPR steps\n"
       << "#NTPL: print out C.O.M motion and total energy fitting every NTPL steps\n"
       << "#NTPP: =1 perform dihedral angle transition monitoring\n"
       << "#     NTPR      NTPL      NTPP\n"
       << setw(10) << gin.print.ntpr
       << setw(10) << gin.print.ntpl
       << setw(10) << gin.print.ntpp
       << "\nEND\n";
  }
  // PRINTOUT (promd, md++)
  if(gin.printout.found){
    os << "PRINTOUT\n"
       << "#NTPR: print out energies, etc. every NTPR steps\n"
       << "#NTPP: =1 perform dihedral angle transition monitoring\n"
       << "#     NTPR      NTPP\n"
       << setw(10) << gin.printout.ntpr
       << setw(10) << gin.printout.ntpp
       << "\nEND\n";
  }
  // WRITE (g96)
  if(gin.write.found){
    os << "WRITE\n"
       << "# NTPW = 0 : binary\n"
       << "# NTPW = 1 : formatted\n"
       << "# NTWSE = configuration selection parameter\n"
       << "# =0: write normal trajectory\n"
       << "# >0: chose min energy for writing configurations\n"
       << "#     NTWX     NTWSE      NTWV      NTWE      NTWG      ";
    if (gin.write.ntba != -1)
      os << "NTBA      NTPW\n";
    else os << "NTPW\n";
    
    os << setw(10) << gin.write.ntwx 
       << setw(10) << gin.write.ntwse
       << setw(10) << gin.write.ntwv 
       << setw(10) << gin.write.ntwe 
       << setw(10) << gin.write.ntwg;
    if (gin.write.ntba != -1)
      os << setw(10) << gin.write.ntba;

    os << setw(10) << gin.write.ntpw
       << "\nEND\n";
  }
  // WRITETRAJ (g96)
  if(gin.writetraj.found){
    os << "WRITETRAJ\n"
       << "# NTPW = 0 : binary\n"
       << "# NTPW = 1 : formatted\n"
       << "# NTWSE = configuration selection parameter\n"
       << "# =0: write normal trajectory\n"
       << "# >0: chose min energy for writing configurations\n"
       << "#    NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB\n"
       << setw(9)  << gin.writetraj.ntwx
       << setw(10) << gin.writetraj.ntwse
       << setw(10) << gin.writetraj.ntwv
       << setw(10) << gin.writetraj.ntwf
       << setw(10) << gin.writetraj.ntwe
       << setw(10) << gin.writetraj.ntwg
       << setw(10) << gin.writetraj.ntwb
       << "\nEND\n";
  }
  // EWARN (md++)
  if(gin.ewarn.found){
    os << "EWARN\n"
       << "#  MAXENER\n"
       << setw(10) << gin.ewarn.maxener
       << "\nEND\n";
  }
  // DEBUG (promd)
  if(gin.debug.found){
    os << "DEBUG\n"
       << "#    NRDEO"
       << setw(10) << gin.debug.routines.size()
       << "\n#  PIIDER(1...NRDEO)  IIIDEO(1...NRDEO)\n";
    for(unsigned int i=0; i<gin.debug.routines.size(); ++i){
      os << setw(25) << gin.debug.routines[i].piider
         << setw(8)  << gin.debug.routines[i].iiideo
         << "\n";
    }
    os << "END\n";
  }
  // SHAKE (g96)
  if(gin.shake.found){
    os << "SHAKE\n"
       << "#      NTC       TOL\n"
       << setw(10) << gin.shake.ntc
       << setw(10) << gin.shake.tol
       << "\nEND\n";
  }
  // GEOMCONSTRAINTS (promd)
  if(gin.geomconstraint.found){
    os << "GEOMCONSTRAINT\n"
       << "#    NTCPH     NTPCH      NTCS    SHKTOL\n"
       << setw(10) << gin.geomconstraint.ntcph
	   << setw(10) << gin.geomconstraint.ntcpn
       << setw(10) << gin.geomconstraint.ntcs
       << setw(10) << gin.geomconstraint.shktol
       << "\nEND\n";
  }
  // CONSTRAINT (md++)
  if(gin.constraint.found){
    os << "CONSTRAINT\n"
       << "# NTC\n"
       << setw(5) << gin.constraint.ntc
       << "\n";
    if(gin.constraint.ntcp == "flexshake"){
      os << "#      NTCP   NTCP0(1 ... 3)\n"
         << setw(11) << gin.constraint.ntcp
         << setw(10) << gin.constraint.ntcp0[0]
         << setw(10) << gin.constraint.ntcp0[1]
         << setw(10) << gin.constraint.ntcp0[2];
    }
    else{
      os << "#      NTCP  NTCP0(1)\n"
         << setw(11) << gin.constraint.ntcp
         << setw(10) << gin.constraint.ntcp0[0];
    }
    os << "\n";
    if(gin.constraint.ntcs == "flexshake"){
      os << "#      NTCS   NTCS0(1 ... 3)\n"
         << setw(11) << gin.constraint.ntcs
         << setw(10) << gin.constraint.ntcs0[0]
         << setw(10) << gin.constraint.ntcs0[1]
         << setw(10) << gin.constraint.ntcs0[2];
    }
    else{
      os << "#      NTCS  NTCS0(1)\n"
         << setw(11) << gin.constraint.ntcs
         << setw(10) << gin.constraint.ntcs0[0];
    }
    os << "\nEND\n";
  }
  // FORCE (promd, md++, g96)
  if(gin.force.found){
    os << "FORCE\n"
       << "#      NTF array\n"
       << "# bonds    angles   imp.     dihe     charge nonbonded\n"
       << "# H        H        H        H\n"
       << setw(3) << gin.force.ntf[0] << setw(3) << gin.force.ntf[1] 
       << setw(6) << gin.force.ntf[2] << setw(3) << gin.force.ntf[3]
       << setw(6) << gin.force.ntf[4] << setw(3) << gin.force.ntf[5]
       << setw(6) << gin.force.ntf[6] << setw(3) << gin.force.ntf[7]
       << setw(6) << gin.force.ntf[8] << setw(3) << gin.force.ntf[9]
       << "\n# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)\n"
       << setw(6) << gin.force.nre.size() << "\n";
	int countnre=0;
    for(unsigned int i=0; i< gin.force.nre.size(); i++){
      os << setw(9) << gin.force.nre[i];
	  countnre++;
	  if(countnre%8==0) os << endl;
    }
    os << "\nEND\n";
  }
  // COVALENTFORM (promd, md++)
  if(gin.covalentform.found){
    os << "COVALENTFORM\n"
       << "#    NTBBH    NTBAH     NTBDN\n"
       << setw(10) << gin.covalentform.ntbbh
       << setw(10) << gin.covalentform.ntbah
       << setw(10) << gin.covalentform.ntbdn
       << "\nEND\n";
  }
  // PLIST03 (old md++ format)
  if(gin.plist03.found){
    os << "PLIST03\n"
       << "#    ALG        NSNB     RCUTP     RCUTL    SIZE      CUTOFF\n";
    if (gin.plist03.grid){
      os << setw(10) << "grid";
    }
    else
      os << setw(10) << "standard";

    os << setw(10) << gin.plist03.nsnb
       << setw(10) << gin.plist03.rcutp
       << setw(10) << gin.plist03.rcutl;
    if (gin.plist03.grda){
      os << setw(10) << "auto";
    }
    else
      os << setw(10) << gin.plist03.grds;

    if (gin.plist03.chargegroup){
      os << setw(15) << "chargegroup";
    }
    else
      os << setw(15) << "atomic";
    os << "\nEND\n";
  }
  // PLIST (g96)
  if(gin.plist.found){
    os << "PLIST\n"
       << "#     NTNB      NSNB     RCUTP     RCUTL\n"
       << setw(10) << gin.plist.ntnb
       << setw(10) << gin.plist.nsnb
       << setw(10) << gin.plist.rcutp
       << setw(10) << gin.plist.rcutl
       << "\nEND\n";
  }
  // NEIGHBOURLIST (promd)
  if(gin.neighbourlist.found){
    os << "NEIGHBOURLIST\n"
       << "# NPLSR  NMESR  NTRSR  RLSTP  RCUTP\n"
       << setw(7) << gin.neighbourlist.nplsr
       << setw(7) << gin.neighbourlist.nmesr
       << setw(7) << gin.neighbourlist.ntrsr
       << setw(7) << gin.neighbourlist.rlstp
       << setw(7) << gin.neighbourlist.rcutp
       << "\n"
       << "# NPLLR  NMELR  NTRLR  RLSTS  RCUTS  RLSTL  RCUTL\n"
       << setw(7) << gin.neighbourlist.npllr
       << setw(7) << gin.neighbourlist.nmelr
       << setw(7) << gin.neighbourlist.ntrlr
       << setw(7) << gin.neighbourlist.rlsts
       << setw(7) << gin.neighbourlist.rcuts
       << setw(7) << gin.neighbourlist.rlstl
       << setw(7) << gin.neighbourlist.rcutl
       << "\n"
       << "# NPLGX  NPLGY  NPLGZ\n"
       << setw(7) << gin.neighbourlist.nplgx
       << setw(7) << gin.neighbourlist.nplgy
       << setw(7) << gin.neighbourlist.nplgz
       << "\nEND\n";
  }
  // PAIRLIST (md++)
  if(gin.pairlist.found){
    os << "PAIRLIST\n"
       << "# algorithm    NSNB   RCUTP   RCUTL    SIZE    TYPE\n"
       << setw(11) << gin.pairlist.algorithm
       << setw(8) << gin.pairlist.nsnb
       << setw(8) << gin.pairlist.rcutp
       << setw(8) << gin.pairlist.rcutl
       << setw(8) << gin.pairlist.size
       << setw(8) << gin.pairlist.type
       << "\nEND\n";
  }
  // NONBONDED (promd)
  if(gin.nonbonded.found){
    os << "NONBONDED\n"
       << "# NLRELE\n"
       << setw(10) << gin.nonbonded.nlrele
       << "\n"
       << "#  APPAK    RCRF   EPSRF\n"
       << setw(10) << gin.nonbonded.appak
       << setw(10) << gin.nonbonded.rcrf
       << setw(10) << gin.nonbonded.epsrf
       << "\n"
       << "# NSHAPE  ASHAPE  NA2CLC   TOLA2   EPSLS\n"
       << setw(10) << gin.nonbonded.nshape
       << setw(10) << gin.nonbonded.ashape
       << setw(10) << gin.nonbonded.na2clc
       << setw(10) << gin.nonbonded.tola2
       << setw(10) << gin.nonbonded.epsls
       << "\n"
       << "#    NKX     NKY     NKZ   KCUT\n"
       << setw(10) << gin.nonbonded.nkx
       << setw(10) << gin.nonbonded.nky
       << setw(10) << gin.nonbonded.nkz
       << setw(10) << gin.nonbonded.kcut
       << "\n"
       << "#    NGX     NGY     NGZ  NASORD  NFDORD  NALIAS  NSPORD\n"
       << setw(10) << gin.nonbonded.ngx
       << setw(10) << gin.nonbonded.ngy
       << setw(10) << gin.nonbonded.ngz
       << setw(10) << gin.nonbonded.nasord
       << setw(10) << gin.nonbonded.nfdord
       << setw(10) << gin.nonbonded.nalias
       << setw(10) << gin.nonbonded.nspord
       << "\n"
       << "# NQEVAL  FACCUR  NRDGRD  NWRGRD\n"
       << setw(10) << gin.nonbonded.nqeval
       << setw(10) << gin.nonbonded.faccur
       << setw(10) << gin.nonbonded.nrdgrd
       << setw(10) << gin.nonbonded.nwrgrd
       << "\n"
       << "#  NLRLJ  SLVDNS\n"
       << setw(10) << gin.nonbonded.nlrlj
       << setw(10) << gin.nonbonded.slvdns
       << "\nEND\n";
  }
  // LONGRANGE (md++, g96)
  if(gin.longrange.found){
    os << "LONGRANGE\n"
       << "# EPSRF     APPAK      RCRF\n"
       << setw(7) << gin.longrange.epsrf
       << setw(10) << gin.longrange.appak
       << setw(10) << gin.longrange.rcrf
       << "\nEND\n";
  }
  // CGRAIN (md++)
  if(gin.cgrain.found){
    os << "CGRAIN\n"
       << "#  NTCGRAN       EPS\n"
       << setw(10) << gin.cgrain.ntcgran
       << setw(10) << gin.cgrain.eps
       << "\nEND\n";
  }
  // POSREST (g96)
  if(gin.posrest.found){
    os << "POSREST\n"
       << "#     values for NTR\n"
       << "#     0: no position re(con)straining\n"
       << "#     1: use CHO\n"
       << "#     2: use CHO/ ATOMIC B-FACTORS\n"
       << "#     3: position constraining\n"
       << "#      NTR       CHO   NRDRX\n"
       << setw(10) << gin.posrest.ntr
       << setw(10) << gin.posrest.cho
       << setw(10) << gin.posrest.nrdrx
       << "\nEND\n";
  }
  // POSITIONRES (promd, md++)
  if(gin.positionres.found){
    os << "POSITIONRES\n"
       << "#     values for NTPOR\n"
       << "#     0: no position re(con)straining\n"
       << "#     1: use CPOR\n"
       << "#     2: use CPOR/ ATOMIC B-FACTORS\n"
       << "#     3: position constraining\n"
       << "#    NTPOR    NTPORB  NTPORS      CPOR\n"
       << setw(10) << gin.positionres.ntpor
       << setw(10) << gin.positionres.ntporb
       << setw(10) << gin.positionres.ntpors
       << setw(10) << gin.positionres.cpor
       << "\nEND\n";
  }
  // DISTREST (g96)
  if(gin.distrest.found){
    os << "DISTREST\n"
       << "#     NTDR\n"
       << "#     0 : no distance restraining\n"
       << "#     -1,1 : use CDIS\n"
       << "#     -2,2:  use W0*CDIS\n"
       << "#     NTDR < 0 : time averaging\n"
       << "#     NTDR > 0 : no time averaging\n"
       << "# NRDDR = 1: read in time averaged distances (for continuation run)\n"
       << "# NRDDR = 0: don't read them in, recalc from scratch\n"
       << "#     NTDR      CDIS       DR0     TAUDR     NRDDR\n"
       << setw(10) << gin.distrest.ntdr
       << setw(10) << gin.distrest.cdis
       << setw(10) << gin.distrest.dr0
       << setw(10) << gin.distrest.taudr
       << setw(10) << gin.distrest.nrddr
       << "\nEND\n";
  }
  // DISTANCERES (promd, md++)
  if(gin.distanceres.found){
    os << "DISTANCERES\n"
       << "#     NTDIR\n"
       << "#     0 : no distance restraining\n"
       << "#     -1,1 : use CDIS\n"
       << "#     -2,2:  use W0*CDIS\n"
       << "#     NTDIR < 0 : time averaging\n"
       << "#     NTDIR > 0 : no time averaging\n"
       << "# NTDIRA = 1: read in time averaged distances (for continuation run)\n"
       << "# NTDIRA = 0: don't read them in, recalc from scratch\n"
       << "#    NTDIR    NTDIRA      CDIR      DIR0    TAUDIR\n"
       << setw(10) << gin.distanceres.ntdir
       << setw(10) << gin.distanceres.ntdira
       << setw(10) << gin.distanceres.cdir
       << setw(10) << gin.distanceres.dir0
       << setw(10) << gin.distanceres.taudir
       << "\nEND\n";
  }
  // DIHEREST (g96)
  if(gin.diherest.found){
    os << "DIHEREST\n"
       << "#    NTDLR      CDLR\n"
       << setw(10) << gin.diherest.ntdlr
       << setw(10) << gin.diherest.cdlr
       << "\nEND\n";
  }
  // DIHEDRALRES (promd,md++)
  if(gin.dihedralres.found){
    os << "DIHEDRALRES\n"
       << "#          NTDLR      CDLR    PHILIN\n"
       << setw(16) << gin.dihedralres.ntdlr
       << setw(10) << gin.dihedralres.cdlr
       << setw(10) << gin.dihedralres.philin
       << "\nEND\n";
  }
  // J-VAL (g96)
  if(gin.jval.found){
    os << "J-VAL\n"
       << "#     NTJR     NTJRH       CJR     TAUJR     NRDJR\n"
       << setw(10) << gin.jval.ntjr
       << setw(10) << gin.jval.ntjrh
       << setw(10) << gin.jval.cjr
       << setw(10) << gin.jval.taujr
       << setw(10) << gin.jval.nrdjr
       << "\nEND\n";
  }
  // JVALUERES (promd, md++)
  if(gin.jvalueres.found){
    os << "JVALUERES\n"
       << "#        NTJVR  NTJVRA    CJVR  TAUJVR\n"
       << setw(16) << gin.jvalueres.ntjvr
       << setw(10) << gin.jvalueres.ntjvra
       << setw(10) << gin.jvalueres.cjvr
       << setw(10) << gin.jvalueres.taujvr
       << "\n"
       << "#     LE   NGRID    DELTA\n"
       << setw(10) << gin.jvalueres.le
       << setw(10) << gin.jvalueres.ngrid
       << setw(10) << gin.jvalueres.delta
       << "\nEND\n";
  }
  // LOCALELEVATION (g96)
  if(gin.localelevation.found){
    os << "LOCALELEVATION\n"
       << "#     NTLE      CWLE     NRDLE\n"
       << setw(10) << gin.localelevation.ntle
       << setw(10) << gin.localelevation.cwle
       << setw(10) << gin.localelevation.nrdle
       << "\nEND\n";
  }
  // LOCALELEV (promd)
  if(gin.localelev.found){
    os << "LOCALELEV\n"
       << "#    NTLES    NTLESA      CLES\n"
       << setw(10) << gin.localelev.ntles
       << setw(10) << gin.localelev.ntlesa
       << setw(10) << gin.localelev.cles
       << "\nEND\n";
  }
  // ROTTRANS (md++)
  if(gin.rottrans.found){
    os << "ROTTRANS\n"
       << "#      RTC   RTCLAST\n"
       << setw(10) << gin.rottrans.rtc
       << setw(10) << gin.rottrans.rtclast
       << "\nEND\n";
  }
  // PERTURB (g96)
  if(gin.perturb.found){
    os << "PERTURB\n"
       << "# NTG: 0 no perturbation is applied\n"
       << "#    : 1 calculate dV/dRLAM perturbation\n"
       << "#    : 2 calculate dV/dRMU perturbation\n"
       << "#    : 3 calculate both derivatives\n"
       << "#      NTG     NRDGL     RLAM     DLAMT        RMU   DMUT\n"
       << setw(10) << gin.perturb.ntg
       << setw(10) << gin.perturb.nrdgl
       << setw(9) << gin.perturb.rlam
       << setw(10) << gin.perturb.dlamt
       << setw(11) << gin.perturb.rmu
       << setw(7) << gin.perturb.dmut
       << "\n#   ALPHLJ     ALPHC     NLAM       MMU\n"
       << setw(10) << gin.perturb.alphlj
       << setw(10) << gin.perturb.alphc
       << setw(9) << gin.perturb.nlam
       << setw(10) << gin.perturb.mmu
       << "\nEND\n";
  }
  // PERTURBATION (promd, md++)
  if(gin.perturbation.found){
    os << "PERTURBATION\n"
       << "# NTG: 0 no perturbation is applied\n"
       << "#    : 1 calculate dV/dRLAM perturbation\n"
       << "#      NTG     NRDGL     RLAM     DLAMT\n"
       << setw(10) << gin.perturbation.ntg
       << setw(10) << gin.perturbation.nrdgl
       << setw(9) << gin.perturbation.rlam
       << setw(10) << gin.perturbation.dlamt
       << "\n#   ALPHLJ     ALPHC     NLAM\n"
       << setw(10) << gin.perturbation.alphlj
       << setw(10) << gin.perturbation.alphc
       << setw(9) << gin.perturbation.nlam
       << "\n#   NSCALE\n"
       << setw(10) << gin.perturbation.nscale
       << "\nEND\n";
  }
  // LAMBDAS (md++)
  if(gin.lambdas.found){
    os << "LAMBDAS\n"
       << "#       NTIL\n"
       << setw(13) << gin.lambdas.lambints.size()
       << "\n# NTLI(1...)  NILG1  NILG2    ALI    BLI    CLI    DLI    ELI\n";
    for(unsigned int i=0; i<gin.lambdas.lambints.size(); ++i){
       os << setw(13) << gin.lambdas.lambints[i].ntli
          << setw(7)  << gin.lambdas.lambints[i].nilg1
          << setw(7)  << gin.lambdas.lambints[i].nilg2
          << setw(7)  << gin.lambdas.lambints[i].ali
          << setw(7)  << gin.lambdas.lambints[i].bli
          << setw(7)  << gin.lambdas.lambints[i].cli
          << setw(7)  << gin.lambdas.lambints[i].dli
          << setw(7)  << gin.lambdas.lambints[i].eli
          << "\n";
    }
    os << "END\n";
  }
  // PERTURB03 (old md++ format)
  if(gin.perturb03.found){
    os << "PERTURB03\n"
       << "# NTG: 0 no perturbation is applied\n"
       << "#    : 1 calculate dV/dRLAM perturbation\n"
       << "#      NTG     RLAM     DLAMT    NLAM   a(LJ)   a(CRF)     SCALING\n"
       << setw(10) << gin.perturb03.ntg
       << setw(9) << gin.perturb03.rlam
       << setw(10) << gin.perturb03.dlamt
       << setw(9) << gin.perturb03.nlam
       << setw(10) << gin.perturb03.alphlj
       << setw(10) << gin.perturb03.alphc
       << setw(10) << gin.perturb03.scaling
       << "\nEND\n";
  }
  // UMBRELLA (promd)
  if(gin.umbrella.found){
    os << "UMBRELLA\n"
       << "#  NTUS  USTC1  USTC2 USREF1 USREF2\n"
       << setw(7) << gin.umbrella.ntus
       << setw(7) << gin.umbrella.uscst1
       << setw(7) << gin.umbrella.uscst2
       << setw(7) << gin.umbrella.usref1
       << setw(7) << gin.umbrella.usref2
       << "\nEND\n";
  }
  // PERSCALE (md++)
  if(gin.perscale.found){
    os << "PERSCALE\n"
       << "#  RESTYPE\n"
       << setw(10) << gin.perscale.restype
       << "\n#     KDIH       KJ         T      DIFF     RATIO      READ\n"
       << setw(10) << gin.perscale.kdih
       << setw(10) << gin.perscale.kj
       << setw(10) << gin.perscale.t
       << setw(10) << gin.perscale.diff
       << setw(10) << gin.perscale.ratio
       << setw(10) << gin.perscale.read
       << "\nEND\n";
  }
  // REPLICA (md++)
  if(gin.replica.found){
    os << "REPLICA\n"
       << "#     NRET\n"
       << setw(10) << gin.replica.ret.size()
       << "\n#  RET(1 ... NRET)\n";
    for(unsigned int i=0; i<gin.replica.ret.size(); ++i){
      os << setw(10) << gin.replica.ret[i];
    }
    os << "\n# LRESCALE\n"
       << gin.replica.lrescale
       << "\n#   NRELAM\n"
       << gin.replica.relam.size()
       << "\n#  RELAM(1 ... NRELAM)\n";
    for(unsigned int i=0; i<gin.replica.relam.size(); ++i){
      os << setw(10) << gin.replica.relam[i];
    }
    os << "\n#   RETS(1 ... NRELAM)\n";
    for(unsigned int i=0; i<gin.replica.rets.size(); ++i){
      os << setw(10) << gin.replica.rets[i];
    }
    os << "\n# NRETRIAL   NREQUIL    NREJOB    NREWRT\n"
       << setw(10) << gin.replica.nretrial
       << setw(10) << gin.replica.nrequil
       << setw(10) << gin.replica.nrejob
       << setw(10) << gin.replica.nrewrt
       << "\nEND\n";
  }
  // INNERLOOP (md++)
  if(gin.innerloop.found){
    os << "INNERLOOP\n"
       << "#     SPEC\n"
       << setw(10) << gin.innerloop.spec
       << "\nEND\n";
  }
  // FOURDIM (g96)
  if(gin.fourdim.found){
    os << "FOURDIM\n"
       << "#  NT4DIM    CW4DA   TEMP4I   NDFMI4    NT4XI    NT4X0\n"
       << setw(9) << gin.fourdim.nt4dim
       << setw(9) << gin.fourdim.cw4da
       << setw(9) << gin.fourdim.temp4i
       << setw(9) << gin.fourdim.ndfmi4
       << setw(9) << gin.fourdim.nt4xi
       << setw(9) << gin.fourdim.nt4xo
       << "\n#    NTT4   TEMP04    TAUT4\n"
       << setw(9) << gin.fourdim.ntt4
       << setw(9) << gin.fourdim.temp04
       << setw(9) << gin.fourdim.taut4
       << "\n#  NTCW4D    CW4DB   TEMP0B\n"
       << setw(9) << gin.fourdim.ntcw4d
       << setw(9) << gin.fourdim.cw4db
       << setw(9) << gin.fourdim.temp0b
       << "\n# NTF4(1)  NTF4(2)  NTF4(3)  NTF4(4)  NTF4(5)  NTF4(6)\n"
       << setw(9) << gin.fourdim.ntf4[0]
       << setw(9) << gin.fourdim.ntf4[1]
       << setw(9) << gin.fourdim.ntf4[2]
       << setw(9) << gin.fourdim.ntf4[3]
       << setw(9) << gin.fourdim.ntf4[4]
       << setw(9) << gin.fourdim.ntf4[5]
       << "\nEND\n";
  }
  // GROMOS96COMPAT (promd)
  if(gin.gromos96compat.found){
    os << "GROMOS96COMPAT\n"
       << "#   NTNB96    NTR96     NTP96     NTG96\n"
	   << setw(10) << gin.gromos96compat.ntnb96
       << setw(10) << gin.gromos96compat.ntr96
       << setw(10) << gin.gromos96compat.ntp96
       << setw(10) << gin.gromos96compat.ntg96
       << "\nEND\n";
  }
  // PATHINT (promd)
  if(gin.pathint.found){
    os << "PATHINT\n"
       << "#  NTPI\n"
       << setw(7) << gin.pathint.ntpi
       << "\nEND\n";
  }
  // INTEGRATE (md++)
  if(gin.integrate.found){
    os << "INTEGRATE\n"
       << "#   NINT\n"
       << setw(8) << gin.integrate.nint
       << "\nEND\n";
  }
  // RANDOMNUMBERS (md++)
  if(gin.randomnumbers.found){
    os << "RANDOMNUMBERS\n"
      << "#   NTRNG   NTGSL\n"
      << setw(8) << gin.randomnumbers.ntrng
	  << setw(8) << gin.randomnumbers.ntgsl
      << "\nEND\n";
  }
  // Unknown blocks
  for(unsigned int i=0; i< gin.unknown.size(); i++){
    os << gin.unknown[i].name << "\n"
       << gin.unknown[i].content
       << "END\n";
  }
    
  return os;
  
}
