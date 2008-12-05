// mk_script.h

// We define two global variables. It makes the printing and bookkeeping 
// of the errors and warnings a lot cleaner.
int numWarnings=0;
int numErrors=0;

enum filetype{unknownfile, inputfile, topofile, coordfile, refposfile, 
	      posresspecfile, disresfile, pttopofile, dihresfile, jvaluefile,
	      ledihfile, outputfile, outtrxfile, outtrvfile, outtrffile, 
	      outtrefile, outtrgfile, 
	      scriptfile, outbaefile, outbagfile};
int numFiletypes=20;
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
		       FT("outtrf", outtrffile),
		       FT("outtre", outtrefile),
		       FT("outtrg", outtrgfile),
		       FT("outbae", outbaefile),
		       FT("outbag", outbagfile),
		       FT("script", scriptfile)
};
static map<string,filetype> FILETYPE(filetypes,filetypes+numFiletypes);

enum blocktype { unknown            , barostatblock    , boundcondblock,
		 cgrainblock        , comtransrotblock , consistencycheckblock,
		 constraintblock    , covalentformblock, debugblock,
		 dihedralresblock   , distanceresblock , energyminblock,
		 ewarnblock         , forceblock       , geomconstraintsblock,
		 gromos96compatblock, initialiseblock  , innerloopblock,
		 integrateblock     , jvalueresblock   , lambdasblock,
		 localelevblock     , multibathblock   , multicellblock,
		 neighbourlistblock , nonbondedblock   , overalltransrotblock,
		 pairlistblock      , pathintblock     , perscaleblock,
		 perturbationblock  , polariseblock    , positionresblock,
		 pressurescaleblock , printoutblock    , randomnumbersblock,
		 readtrajblock      , replicablock     , rottransblock,
		 stepblock          , stochdynblock    , systemblock,
		 thermostatblock    , umbrellablock    , virialblock,
		 writetrajblock };

typedef map<string, blocktype>::value_type BT;
int numBlocktypes = 46;
const BT blocktypes[] ={BT("",unknown),
			BT("BAROSTAT",barostatblock),			
			BT("BOUNDCOND",boundcondblock),
			BT("CGRAIN",cgrainblock),
			BT("COMTRANSROT",comtransrotblock),
			BT("CONSISTENCYCHECK",consistencycheckblock),
			BT("CONSTRAINT",constraintblock),
			BT("COVALENTFORM",covalentformblock),
			BT("DEBUG",debugblock),
			BT("DIHEDRALRES",dihedralresblock),
			BT("DISTANCERES",distanceresblock),
			BT("ENERGYMIN",energyminblock),
			BT("EWARN",ewarnblock),
			BT("FORCE",forceblock),
			BT("GEOMCONSTRAINTS",geomconstraintsblock),
			BT("GROMOS96COMPAT",gromos96compatblock),
			BT("INITIALISE",initialiseblock),
			BT("INNERLOOP",innerloopblock),
			BT("INTEGRATE",integrateblock),
			BT("JVALUERES",jvalueresblock),
			BT("LAMBDAS", lambdasblock),
			BT("LOCALELEV",localelevblock),
			BT("MULTIBATH",multibathblock),
			BT("MULTICELL",multicellblock),
			BT("NEIGHBOURLIST",neighbourlistblock),
			BT("NONBONDED",nonbondedblock),
			BT("OVERALLTRANSROT",overalltransrotblock),
			BT("PAIRLIST",pairlistblock),
			BT("PATHINT",pathintblock),
			BT("PERSCALE",perscaleblock),
			BT("PERTURBATION",perturbationblock),
			BT("POLARISE", polariseblock),
			BT("POSITIONRES",positionresblock),
			BT("PRESSURESCALE",pressurescaleblock),
			BT("PRINTOUT",printoutblock),
			BT("RANDOMNUMBERS",randomnumbersblock),
			BT("READTRAJ",readtrajblock),
			BT("REPLICA",replicablock),
			BT("ROTTRANS",rottransblock),
			BT("STEP",stepblock),
			BT("STOCHDYN",stochdynblock),
      			BT("SYSTEM",systemblock),
			BT("THERMOSTAT",thermostatblock),
			BT("UMBRELLA",umbrellablock),
			BT("VIRIAL",virialblock),
			BT("WRITETRAJ",writetrajblock),
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

class iboundcond{
public:
  int found, ntb, ndfmin;
  iboundcond(){found=0;}
};

class icgrain{
public:
  int found, ntcgran;
  double eps;
  icgrain(){found=0;}
};

class icomtransrot{
public:
  int found, nscm;
  icomtransrot(){found=0;}
};
 
class iconsistencycheck{
public:
  int found, ntchk, ntckf, ntckv, ntckt;
  int ntcke, ntckr, ntckl;
  double fdckf, fdckv, fdckl;
  vector<int> nckf;
  iconsistencycheck(){found=0;}
};

class iconstraint{
public:
  int found, ntc, ntcp, ntcs;
  double ntcp0[3], ntcs0[3]; 
  iconstraint(){found=0;}
};

class icovalentform{
public:
  int found, ntbbh, ntbah, ntbdn;
  icovalentform(){found=0;}
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

class idihedralres{
public:
  int found, ntdlr;
  double cdlr, philin;
  idihedralres(){found=0;}
};

class idistanceres{
public:
  int found, ntdir, ntdira;
  double cdir, dir0, taudir;
  idistanceres(){found=0;}
};

class ienergymin{
public:
  int found, ntem, ncyc, nmin;
  double dele, dx0, dxm, flim;
  ienergymin(){found=0;}
};

class iewarn{
public:
  int found;
  double maxener;
  iewarn(){found=0;}
};

class iforce{
public:
  int found, ntf[10];
  vector<int> nre;
  iforce(){found=0;}
};

class igeomconstraints{
public:
  int found, ntcph, ntcpn, ntcs;
  double shktol;
  igeomconstraints(){found=0;}
}; 

class igromos96compat{
public:
  int found, ntnb96, ntr96, ntp96, ntg96;
  igromos96compat(){found=0;}
};

class iinitialise{
public:
  int found, ntivel, ntishk, ntinht, ntinhb;
  int ntishi, ntirtc, nticom, ntisti;
  double ig, tempi;
  iinitialise(){found=0;}
};

class iinnerloop{
public:
  int found, ntilm, ntils;
  iinnerloop(){found=0;}
};

class iintegrate{
public:
  int found, nint;
  iintegrate(){found=0;}
};

class ijvalueres{
public:
  int found, ntjvr, ntjvra, le, ngrid;
  double cjvr, taujvr, delta;
  ijvalueres(){found=0;}
};

class ilambdas{
public:
  int found, ntil;
  class lambint{
  public:
    int ntli, nilg1, nilg2;
    double ali, bli, cli, dli, eli;
  };
  vector<lambint> lambints;
  ilambdas(){found=0;}
};

class ilocalelev{
public:
  int found, ntles, ntlefr, ntlefu, nlegrd, ntlesa;
  double cles, wles, rles;
  ilocalelev(){found=0;}
};

class imultibath{
public:
  int found, algorithm, num, nbaths, dofset;
  vector<double> temp0, tau;
  vector<int> last, combath, irbath;
  imultibath(){found=0; num=-1;}
};

class imulticell{
public:
  int found, ntm, ncella, ncellb, ncellc;
  double tolpx, tolpv, tolpf, tolpfw;
  imulticell(){found=0;}
};

class ineighbourlist{
public:
  int found, nmprpl, nuprpl,nmtwpl, nutwpl, nuirin, nusrin, nmtwin, ncgcen;
  double rcprpl, grprpl, rstwpl, rltwpl, rctwin;
  ineighbourlist(){found=0;}
};

class inonbonded{
public:
  int found, nlrele, nshape, na2clc, nkx, nky, nkz, ngx, ngy, ngz;
  int nasord, nfdord, nalias, nqeval, nrdgrd, nwrgrd, nlrlj;
  double appak, rcrf, epsrf, ashape, tola2, epsls, kcut, nspord, faccur, slvdns;
  inonbonded(){found=0;}
};

class ioveralltransrot{
public:
  int found, ncmtr, ncmro;
  double cmamx, cmamy, cmamz;
  ioveralltransrot(){found=0;}
};

class ipairlist{
public:
  int found, algorithm, nsnb, type;
  double rcutp, rcutl;
  string size;
  ipairlist(){found=0;}
};

class ipathint{
public:
  int found, ntpi;
  ipathint(){found=0;}
};

class iperscale{
public:
  int found, t, read;
  double kdih, kj, diff, ratio;
  string restype;
  iperscale(){found=0;}
};

class iperturbation{
public:
  int found, ntg, nrdgl, nlam, nscale;
  double rlam, dlamt, alphlj, alphc;
  iperturbation(){found=0;}
};

class ipolarise{
public:
  int found, cos, efield, damp, write;
  double minfield;
  ipolarise(){found=0;}
};

class ipositionres{
public:
  int found, ntpor, ntporb, ntpors;
  double cpor;
  ipositionres(){found=0;}
};

class ipressurescale{
public:
  int found, couple, scale, virial;
  double comp, taup, pres0[3][3];
  ipressurescale(){found=0;}
};

class iprintout{
public:
  int found, ntpr, ntpp;
  iprintout(){found=0;}
};

class irandomnumbers{
public:
  int found, ntrng, ntgsl;
  irandomnumbers(){found=0;}
};

class ireadtraj{
public:
  int found, ntrd, ntrn, ntrb, ntshk;
  ireadtraj(){found=0;}
};

class ireplica{
public:
  int found, lrescale, nretrial, nrequil, nrejob, nrewrt;
  vector<double> ret, relam, rets; 
  ireplica(){found=0; nrewrt=0;}
};

class irottrans{
public:
  int found, rtc, rtclast;
  irottrans(){found=0;}
};

class istep{
public:
  int found, nstlim;
  double t,dt;
  istep(){found=0;}
};

class istochdyn{
public:
  int found, ntsd, ntfr, nsfr, nbref;
  double rcutf, cfric, tempsd;
  istochdyn(){found=0;}
};

class isystem{
public:
  int found, npm, nsm;
  isystem(){found=0;}
};

class ithermostat{
public:
  int found, ntt;
  class tbath{
  public:
    int ntbtyp;
    double tembth;
    vector<double> taubth;
  };
  vector<tbath> baths; 
  class dofgroup{
  public:
    int ntscpl, ntstyp, ntscns, ntsgt;
    vector<int> ntsgtg;
  };
  vector<dofgroup> dofgroups;
  ithermostat(){found=0;}
};

class iumbrella{
public:
  int found, ntus;
  double uscst1, uscst2, usref1, usref2;
  iumbrella(){found=0;}
};

class ivirial{
public:
  int found, ntv, ntvg;
  ivirial(){found=0;}
};

class iwritetraj{
public:
  int found, ntwx, ntwse, ntwv, ntwf, ntwe, ntwg, ntwb;
  iwritetraj(){found=0; ntwx=0; ntwse=0; ntwv=0; ntwf=0; ntwe=0; ntwg=0; ntwb=0;}
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

  ibarostat barostat;
  iboundcond boundcond;
  icgrain cgrain;
  icomtransrot comtransrot;
  iconsistencycheck consistencycheck;
  iconstraint constraint;
  icovalentform covalentform;
  idebug debug;
  idihedralres dihedralres;
  idistanceres distanceres;
  ienergymin energymin;
  iewarn ewarn;
  iforce force;
  igeomconstraints geomconstraints;
  igromos96compat gromos96compat;
  iinitialise initialise;
  iinnerloop innerloop;
  iintegrate integrate;
  ijvalueres jvalueres;
  ilambdas lambdas;
  ilocalelev localelev;
  imultibath multibath;
  imulticell multicell;
  ineighbourlist neighbourlist;
  inonbonded nonbonded;
  ioveralltransrot overalltransrot;
  ipairlist pairlist;
  ipathint pathint;
  iperscale perscale;
  iperturbation perturbation;
  ipolarise polarise;
  ipositionres positionres;
  ipressurescale pressurescale;
  iprintout printout;
  irandomnumbers randomnumbers;
  ireadtraj readtraj;
  ireplica replica;
  irottrans rottrans;
  istep step;
  istochdyn stochdyn;
  isystem system;
  ithermostat thermostat;
  iumbrella umbrella;
  ivirial virial;
  iwritetraj writetraj;
  vector<iunknown> unknown;
}; 

class fileInfo{
public:
  gcore::Box box;
  vector<string> blocks;
  vector<int> blockslength;
};

//INSTREAM
istringstream &operator>>(istringstream &is,ibarostat &s){
  string e;
  int dum, npbth;
  double taup;
  s.found=1;
  is >> s.ntp >> s.npvar >> npbth >> s.comp;
  if(s.ntp!=0 && s.ntp!=1 && s.ntp!=2 && s.ntp!=3){
    std::stringstream ss; ss << s.ntp; 
    printIO("BAROSTAT", "NTP", ss.str(), "0,1,2,3");
  }
  if(s.npvar < 0){
    std::stringstream ss;
    ss << s.npvar;
    printIO("BAROSTAT", "NPVAR", ss.str(), ">= 0");
  }
  if(npbth < 0){
    std::stringstream ss;
    ss << npbth;
    printIO("BAROSTAT", "NPBTH", ss.str(), ">= 0");
  }
  if(s.comp <= 0.0){
    std::stringstream ss;
    ss << s.comp;
    printIO("BAROSTAT", "COMP", ss.str(), "> 0.0");
  }
      
  for(int i=0;i<npbth;++i){
    class ibarostat::pbath p;
    is >> dum >> p.prsbth;
    if(p.prsbth < 0.0){
      std::stringstream si;
      si << "PRSBTH[" << dum << "]";
      std::stringstream ss;
      ss << p.prsbth;
      printIO("BAROSTAT", si.str(), ss.str(), ">= 0.0");
    }
    for(int j=0;j<s.npvar; ++j){
      is >> taup;
      p.taubba.push_back(taup);
      if(taup <=0){
	std::stringstream si;
	si << "TAUBBA[" << dum << "," << j+1 << "]";
	std::stringstream ss;
	ss << taup;
	printIO("BAROSTAT", si.str(), ss.str(), "> 0.0");
      }
    }
    s.pbaths.push_back(p);
  }
  for(int i=0; i< 6; i++){
    is >> s.npcpl[i];
    if(s.npcpl[i] !=0 && s.npcpl[i]!=1){
      std::stringstream si;
      si << "NPCPL[" << i+1 << "]";
      std::stringstream ss;
      ss << s.npcpl[i];
      printIO("BAROSTAT", si.str(), ss.str(), "0,1");
    }
  }
  is >> e;

  return is;
}
istringstream &operator>>(istringstream &is, iboundcond &s){
  string e;
  s.found=1;
  is >> s.ntb >> s.ndfmin >> e;
  if(s.ntb!=-1 && s.ntb!=0 && s.ntb!=1 &&s.ntb!=2){
    std::stringstream ss;
    ss << s.ntb;
    printIO("BOUNDCOND", "NTB", ss.str(), "-1,0,1,2");
  }
  if(s.ndfmin < 0){
    std::stringstream ss;
    ss << s.ndfmin;
    printIO("BOUNDCOND", "NDFMIN", ss.str(), ">=0");
  }
  
  return is;
}
istringstream &operator>>(istringstream &is,icgrain &s){
  string e;
  s.found=1;
  is >> s.ntcgran >> s.eps >> e;
  if(s.ntcgran < 0 || s.ntcgran > 2){
    std::stringstream ss;
    ss << s.ntcgran;
    printIO("CGRAIN", "NTCGRAN", ss.str(), "0,1,2");
  }
  if(s.eps < 0.0){
    std::stringstream ss;
    ss << s.eps;
    printIO("CGRAIN", "EPS", ss.str(), ">= 0.0");
  }
  return is;
}
istringstream &operator>>(istringstream &is,icomtransrot &s){
  string e;
  s.found=1;
  is >> s.nscm >> e;
  // NSCM can take any value
  return is;
}
istringstream &operator>>(istringstream &is, iconsistencycheck &s){
  string e;
  s.found=1;
  is >> s.ntchk >> s.ntckf >> s.fdckf >> s.ntckv >> s.fdckv >> s.ntckt
     >> s.ntcke >> s.ntckr >> s.ntckl >> s.fdckl;
  int nackf, nckf;
  is >> nackf;
  if(nackf < 0){
    std::stringstream ss;
    ss << nackf;
    printIO("CONSISTENCYCHECK", "NACKF", ss.str(), ">= 0");
  }
  for(int i=0; i<nackf; i++) {
    is >> nckf; 
    if(nckf<=0){
      std::stringstream si, ss;
      si << "NCKF[" << i+1 << "]";
      ss << nckf;
      printIO("CONSISTENCYCHECK", si.str(), ss.str(), "> 0");
    }
    s.nckf.push_back(nckf);
  }
  is >> e;
  if(s.ntchk!=0 && s.ntchk!=1){
    std::stringstream ss;
    ss << s.ntchk;
    printIO("CONSISTENCYCHECK", "NTCHK", ss.str(),"0,1");
  }
  if(s.ntckf!=0 && s.ntckf!=1){
    std::stringstream ss;
    ss << s.ntckf;
    printIO("CONSISTENCYCHECK", "NTCKF", ss.str(),"0,1");
  }
  if(s.fdckf<=0.0){
    std::stringstream ss;
    ss << s.fdckf;
    printIO("CONSISTENCYCHECK", "NDCKF", ss.str(),"> 0.0");
  }
  if(s.ntckv!=0 && s.ntckv!=1){
    std::stringstream ss;
    ss << s.ntckv;
    printIO("CONSISTENCYCHECK", "NTCKV", ss.str(),"0,1");
  }
  if(s.fdckv<=0.0){
    std::stringstream ss;
    ss << s.fdckv;
    printIO("CONSISTENCYCHECK", "NDCKV", ss.str(),"> 0.0");
  }
  if(s.ntckt!=0 && s.ntckt!=1){
    std::stringstream ss;
    ss << s.ntckt;
    printIO("CONSISTENCYCHECK", "NTCKT", ss.str(),"0,1");
  }
  if(s.ntcke!=0 && s.ntcke!=1){
    std::stringstream ss;
    ss << s.ntcke;
    printIO("CONSISTENCYCHECK", "NTCKE", ss.str(),"0,1");
  }
  if(s.ntckr!=0 && s.ntckr!=1){
    std::stringstream ss;
    ss << s.ntckr;
    printIO("CONSISTENCYCHECK", "NTCKR", ss.str(),"0,1");
  }
  if(s.ntckl!=0 && s.ntckl!=1){
    std::stringstream ss;
    ss << s.ntckl;
    printIO("CONSISTENCYCHECK", "NTCKL", ss.str(),"0,1");
  }
  if(s.fdckl<=0.0){
    std::stringstream ss;
    ss << s.fdckl;
    printIO("CONSISTENCYCHECK", "NDCKF", ss.str(),"> 0.0");
  }


  return is;
}

istringstream &operator>>(istringstream &is,iconstraint &s){
  string e, ntc, cp, cs;
  s.found=1;
  is >> ntc;
  std::transform(ntc.begin(), ntc.end(), ntc.begin(), ::tolower);
  if(ntc=="solvent") s.ntc=1;
  else if(ntc=="hydrogen") s.ntc=2;
  else if(ntc=="all") s.ntc=3;
  else if(ntc=="specified") s.ntc=4;
  else{
    std::stringstream ss(ntc);
    if((!(ss >> s.ntc)) ||
       (s.ntc != 1 && s.ntc != 2 && s.ntc != 3 && s.ntc !=4))
      printIO("CONSTRAINT", "NTC", ntc, 
	      "solvent(1), hydrogen(2), all(3), specified(4)");
  }

  is>> cp;
  std::transform(cp.begin(), cp.end(), cp.begin(), ::tolower);
  if(cp=="shake") s.ntcp=1;
  else if(cp=="lincs") s.ntcp=2;
  else if(cp=="flexshake") s.ntcp=3;
  else {
    std::stringstream ss(cp);
    if((!(ss >> s.ntcp)) || s.ntcp < 0 || s.ntcp > 3)
      printIO("CONSTRAINT", "NTCP", cp, "shake(1), lincs(2), flexshake(3)");
  }
  is >> s.ntcp0[0];
  if(s.ntcp0[0] < 0){
    std::stringstream ss;
    ss << s.ntcp0[0];
    printIO("CONSTRAINT", "NTCP0[0]", ss.str(), " >= 0.0");
  }
  
  if(s.ntcp==3){
    is >> s.ntcp0[1] >> s.ntcp0[2];
    if(s.ntcp0[1] < 0){
      std::stringstream ss;
      ss << s.ntcp0[1];
      printIO("CONSTRAINT", "NTCP0[1]", ss.str(), " >= 0.0");
    }
    if(s.ntcp0[2] < 0){
      std::stringstream ss;
      ss << s.ntcp0[2];
      printIO("CONSTRAINT", "NTCP0[2]", ss.str(), " >= 0.0");
    }
  }
  is >> cs;
  std::transform(cs.begin(), cs.end(), cs.begin(), ::tolower);
  if(cs=="shake") s.ntcs=1;
  else if(cs=="lincs") s.ntcs=2;
  else if(cs=="settle") s.ntcs=3;
  else{
    std::stringstream ss(cs);
    if((!(ss >> s.ntcs)) || s.ntcs < 0 || s.ntcs > 3)
      printIO("CONSTRAINT", "NTCS", cp, "shake(1), lincs(2), settle(3)");
  }
  if(s.ntcs!=3){
    is >> s.ntcs0[0];
    if(s.ntcs0[0] < 0){
      std::stringstream ss;
      ss << s.ntcs0[0];
      printIO("CONSTRAINT", "NTCS0[0]", ss.str(), " >= 0.0");
    }
  }
  
  is >> e;
  return is;
}

istringstream &operator>>(istringstream &is,icovalentform &s){
  string e;
  s.found=1;
  is >> s.ntbbh >> s.ntbah >> s.ntbdn >> e;
  if(s.ntbbh!=0 && s.ntbbh!=1){
    std::stringstream ss;
    ss << s.ntbbh;
    printIO("COVALENTFORM", "NTBBH", ss.str(), "0,1");
  }
  if(s.ntbah!=0 && s.ntbah!=1){
    std::stringstream ss;
    ss << s.ntbah;
    printIO("COVALENTFORM", "NTBAH", ss.str(), "0,1");
  }
  if(s.ntbdn!=0 && s.ntbdn!=1){
    std::stringstream ss;
    ss << s.ntbdn;
    printIO("COVALENTFORM", "NTBDN", ss.str(), "0,1");
  }
  
  return is;
}

istringstream &operator>>(istringstream &is,idebug &s){
  string e;
  int nrd;
  s.found=1;
  is >> nrd;
  if(nrd<0){
    std::stringstream ss;
    ss << nrd;
    printIO("DEBUG", "NRD", ss.str(), ">= 0");
  }
  for(int i=0; i<nrd; ++i){
    class idebug::routine r;
    is >> r.piider;
    is >> r.iiideo;
    if(r.iiideo <=0){
      std::stringstream ss, si;
      ss << r.iiideo;
      si << "IIIDEO[" << i+1 << "]";
      printIO("DEBUG", si.str(), ss.str(), "> 0");
    }
    s.routines.push_back(r);
  }
  is >> e;
  return is;
}

istringstream &operator>>(istringstream &is,idihedralres &s){
  string e;
  s.found=1;
  is >> s.ntdlr >> s.cdlr >> s.philin >> e;
  if(s.ntdlr <0 || s.ntdlr > 3){
    std::stringstream ss;
    ss << s.ntdlr;
    printIO("DIHEDRALRES", "NTDLR", ss.str(), "0,1,2,3");
  }
  if(s.cdlr < 0.0){
    std::stringstream ss;
    ss << s.cdlr;
    printIO("DIHEDRALRES", "CDLR", ss.str(), ">= 0.0");
  }
  
  return is;
}

istringstream &operator>>(istringstream &is,idistanceres &s){
  string e;
  s.found=1;
  is >> s.ntdir >> s.ntdira >> s.cdir >> s.dir0 >> s.taudir >> e;
  if(s.ntdir < -2 || s.ntdir > 3){
    std::stringstream ss;
    ss << s.ntdir;
    printIO("DISTANCERES", "NTDIR", ss.str(), "-2,-1,0,1,2,3");
  }
  if(s.ntdira!=0 && s.ntdira !=1){
    std::stringstream ss;
    ss << s.ntdir;
    printIO("DISTANCERES", "NTDIRA", ss.str(), "0,1");
  }
  if(s.cdir < 0.0){
    std::stringstream ss;
    ss << s.cdir;
    printIO("DISTANCERES", "CDIR", ss.str(), ">= 0.0");
  }
  if(s.dir0 < 0.0){
    std::stringstream ss;
    ss << s.dir0;
    printIO("DISTANCERES", "DIR0", ss.str(), ">= 0.0");
  }
  if(s.taudir < 0.0){
    std::stringstream ss;
    ss << s.taudir;
    printIO("DISTANCERES", "TAUDIR", ss.str(), ">= 0.0");
  }

  return is;
}

istringstream &operator>>(istringstream &is, ienergymin &s){
  string e;
  s.found=1;
  is >> s.ntem >> s.ncyc >> s.dele >> s.dx0 >> s.dxm >> s.nmin >> s.flim >> e;
  if(s.ntem < 0 || s.ntem > 2){
    std::stringstream ss;
    ss << s.ntem;
    printIO("ENERGYMIN", "NTEM", ss.str(), "0,1,2");
  }
  if(s.ncyc <=0){
    std::stringstream ss;
    ss << s.ncyc;
    printIO("ENERGYMIN", "NCYC", ss.str(), "> 0");
  }
  if(s.dele <= 0.0){
    std::stringstream ss;
    ss << s.dele;
    printIO("ENERGYMIN", "DELE", ss.str(), "> 0.0");
  }
  if(s.dx0 <= 0.0){
    std::stringstream ss;
    ss << s.dx0;
    printIO("ENERGYMIN", "DX0", ss.str(), "> 0.0");
  }
  if(s.dxm <= 0.0){
    std::stringstream ss;
    ss << s.dxm;
    printIO("ENERGYMIN", "DXM", ss.str(), "> 0.0");
  }
  if(s.nmin <= 0){
    std::stringstream ss;
    ss << s.nmin;
    printIO("ENERGYMIN", "NMIN", ss.str(), "> 0");
  }
  if(s.flim < 0.0){
    std::stringstream ss;
    ss << s.flim;
    printIO("ENERGYMIN", "FLIM", ss.str(), ">= 0.0");
  }
  
  return is;
}

istringstream &operator>>(istringstream &is,iewarn &s){
  string e;
  s.found=1;
  is >> s.maxener >> e;

  return is;
}

istringstream &operator>>(istringstream &is,iforce &s){
  string e;
  s.found=1;
  int negr, nre;
  for(int i=0; i<10; i++){
    is >> s.ntf[i];
    if(s.ntf[i]!=0 && s.ntf[i]!=1){
      std::stringstream ss, si;
      ss << s.ntf[i];
      si << "NTF[" << i+1 << "]";
      printIO("FORCE", si.str(), ss.str(), "0,1");
    }
  }
  is >> negr;
  if(negr <= 0){
    std::stringstream ss;
    ss << negr;
    printIO("FORCE", "NEGR", ss.str(), "> 0");
  }
  
  for(int i=0; i<negr; i++) {
    is >> nre; s.nre.push_back(nre);
    if(nre<=0){
      std::stringstream ss, si;
      si << "NRE[" << i+1 << "]";
      ss << nre;
      printIO("FORCE", si.str(), ss.str(), "> 0");
    }
  }
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is,igeomconstraints &s){
  string e;
  s.found=1;
  is >> s.ntcph >> s.ntcpn >> s.ntcs >> s.shktol >> e;
  if(s.ntcph!=0 && s.ntcph!=1){
    std::stringstream ss;
    ss << s.ntcph;
    printIO("GEOMCONSTRAINTS", "NTCPH", ss.str(), "0,1");
  }
  if(s.ntcpn!=0 && s.ntcpn!=1){
    std::stringstream ss;
    ss << s.ntcpn;
    printIO("GEOMCONSTRAINTS", "NTCPN", ss.str(), "0,1");
  }
  if(s.ntcs!=0 && s.ntcs!=1){
    std::stringstream ss;
    ss << s.ntcs;
    printIO("GEOMCONSTRAINTS", "NTCS", ss.str(), "0,1");
  }
  if(s.shktol<=0.0){
    std::stringstream ss;
    ss << s.shktol;
    printIO("GEOMCONSTRAINTS", "SHKTOL", ss.str(), "> 0.0");
  }
  
  return is;
}

istringstream &operator>>(istringstream &is,igromos96compat &s){
  string e;
  s.found=1;
  is >> s.ntnb96 >> s.ntr96 >> s.ntp96 >> s.ntg96 >> e;
  if(s.ntnb96!=0 && s.ntnb96!=1){
    std::stringstream ss;
    ss << s.ntnb96;
    printIO("GROMOS96COMPAT", "NTNB96", ss.str(), "0,1");
  }
  if(s.ntr96!=0 && s.ntr96!=1){
    std::stringstream ss;
    ss << s.ntr96;
    printIO("GROMOS96COMPAT", "NTR96", ss.str(), "0,1");
  }
  if(s.ntp96!=0 && s.ntp96!=1){
    std::stringstream ss;
    ss << s.ntp96;
    printIO("GROMOS96COMPAT", "NTP96", ss.str(), "0,1");
  }
  if(s.ntg96!=0 && s.ntg96!=1){
    std::stringstream ss;
    ss << s.ntg96;
    printIO("GROMOS96COMPAT", "NTG96", ss.str(), "0,1");
  }

  return is;
}

istringstream &operator>>(istringstream &is, iinitialise &s){
  string e;
  s.found=1;
  is >> s.ntivel >> s.ntishk >> s.ntinht >> s.ntinhb;
  is >> s.ntishi >> s.ntirtc >> s.nticom;
  is >> s.ntisti >> s.ig >> s.tempi >> e;
  
  if(s.ntivel!=0 && s.ntivel!=1){
    std::stringstream ss;
    ss << s.ntivel;
    printIO("INITIALISE", "NTIVEL", ss.str(), "0,1");
  }
  if(s.ntishk<0 || s.ntivel>3){
    std::stringstream ss;
    ss << s.ntishk;
    printIO("INITIALISE", "NTISHK", ss.str(), "0,1,2,3");
  }
  if(s.ntinht!=0 && s.ntinht!=1){
    std::stringstream ss;
    ss << s.ntinht;
    printIO("INITIALISE", "NTINHT", ss.str(), "0,1");
  }
  if(s.ntinhb!=0 && s.ntinhb!=1){
    std::stringstream ss;
    ss << s.ntinhb;
    printIO("INITIALISE", "NTINHB", ss.str(), "0,1");
  }
  if(s.ntishi!=0 && s.ntishi!=1){
    std::stringstream ss;
    ss << s.ntishi;
    printIO("INITIALISE", "NTISHI", ss.str(), "0,1");
  }
  if(s.ntirtc!=0 && s.ntirtc!=1){
    std::stringstream ss;
    ss << s.ntirtc;
    printIO("INITIALISE", "NTIRTC", ss.str(), "0,1");
  }
  if(s.nticom<0 || s.nticom > 3){
    std::stringstream ss;
    ss << s.nticom;
    printIO("INITIALISE", "NTICOM", ss.str(), "0,1,2,3");
  }
  if(s.ntisti!=0 && s.ntisti!=1){
    std::stringstream ss;
    ss << s.ntisti;
    printIO("INITIALISE", "NTISTI", ss.str(), "0,1");
  }
  if(s.ig <= 0){
    std::stringstream ss;
    ss << s.ig;
    printIO("INITIALISE", "IG", ss.str(), "> 0");
  }
  if(s.tempi < 0.0){
    std::stringstream ss;
    ss << s.tempi;
    printIO("INITIALISE", "TEMPI", ss.str(), ">= 0");
  }

  return is;
}

istringstream &operator>>(istringstream &is,iinnerloop &s){
  string e;
  s.found=1;
  is >> s.ntilm >> s.ntils >> e;
  if(s.ntilm < 0 || s.ntilm > 4){
    std::stringstream ss;
    ss << s.ntilm;
    printIO("INNERLOOP", "NTILM", ss.str(), "0,1,2,3,4");
  }
  if(s.ntils !=0 && s.ntils !=1){
    std::stringstream ss;
    ss << s.ntils;
    printIO("INNERLOOP", "NTILS", ss.str(), "0,1");
  }
  
  return is;
}

istringstream &operator>>(istringstream &is,iintegrate &s){
  string e;
  s.found=1;
  is >> s.nint >> e;
  if(s.nint!=0 && s.nint!=1){
    std::stringstream ss;
    ss << s.nint;
    printIO("INTEGRATE", "NINT", ss.str(), "0,1");
  }
  return is;
}

istringstream &operator>>(istringstream &is,ijvalueres &s){
  string e;
  s.found=1;
  is >> s.ntjvr >> s.ntjvra >> s.cjvr >> s.taujvr >> s.le 
     >> s.ngrid >> s.delta >> e;
  if(s.ntjvr < -3 || s.ntjvr > 2){
    std::stringstream ss;
    ss << s.ntjvr;
    printIO("JVALUERES", "NTJVR", ss.str(), "-3,-2,-1,0,1,2");
  }
  if(s.ntjvra !=0 && s.ntjvra !=1){
    std::stringstream ss;
    ss << s.ntjvra;
    printIO("JVALUERES", "NTJVRA", ss.str(), "0,1");
  }
  if(s.cjvr < 0.0){
    std::stringstream ss;
    ss << s.cjvr;
    printIO("JVALUERES", "CJVR", ss.str(), ">= 0.0");
  }
  if(s.taujvr < 0.0){
    std::stringstream ss;
    ss << s.taujvr;
    printIO("JVALUERES", "TAUJVR", ss.str(), ">= 0.0");
  }
  if(s.le !=0 && s.le !=1){
    std::stringstream ss;
    ss << s.le;
    printIO("JVALUERES", "LE", ss.str(), "0,1");
  }
  if(s.ngrid <= 0){
    std::stringstream ss;
    ss << s.ngrid;
    printIO("JVALUERES", "NGRID", ss.str(), "> 0");
  }
  if(s.delta < 0.0){
    std::stringstream ss;
    ss << s.delta;
    printIO("JVALUERES", "DELTA", ss.str(), ">= 0.0");
  }

  return is;
}

istringstream &operator>>(istringstream &is,ilambdas &s){
  string e, dum;
  s.found=1;
  is >> dum;
  std::transform(dum.begin(), dum.end(), dum.begin(), ::tolower);
  if(dum=="on") s.ntil=1;
  else if(dum=="off") s.ntil=0;
  else{
    std::stringstream ss(dum);
    if(!(ss>>s.ntil) || (s.ntil !=0 && s.ntil !=1)){
      std::stringstream ss;
      ss << s.ntil;
      printIO("LAMBDAS", "NTIL", ss.str(), "off(0),on(1)");
    }
  }
  is >> dum;
  int i=0;
  while(dum!="END"){
    i++;
    std::transform(dum.begin(), dum.end(), dum.begin(), ::tolower);
    class ilambdas::lambint l;
    if(dum=="bond") l.ntli=1;
    else if(dum=="angle") l.ntli=2;
    else if(dum=="dihedral") l.ntli=3;
    else if(dum=="improper") l.ntli=4;
    else if(dum=="vdw") l.ntli=5;
    else if(dum=="vdw_soft")  l.ntli=6;
    else if(dum=="crf") l.ntli=7;
    else if(dum=="crf_soft") l.ntli=8;
    else if(dum=="distanceres") l.ntli=9;
    else if(dum=="dihedralres") l.ntli=10;
    else if(dum=="mass") l.ntli=11;
    else{
      std::stringstream ss(dum);
      if((!(ss >> l.ntli)) || l.ntli < 1 || l.ntli > 11){
	std::stringstream si;
	si << "NTLI[" << i << "]";
	printIO("LAMBDAS", si.str(), dum, 
		"bond(1), angle(2), dihedral(3), improper(4), vdw(5), "
		"vdw_soft(6), crf(7), crf_soft(8), distanceres(9), "
		"dihedralres(10), mass(11)");
      }
    }
    is >> l.nilg1 >> l.nilg2 >> l.ali >> l.bli >> l.cli >> l.dli >> l.eli;
    if(l.nilg1 <= 0){
      std::stringstream ss, si;
      ss << l.nilg1;
      si << "NILG1[" << i << "]";
      printIO("LAMBDAS", si.str(), ss.str(), "> 0");
    }
    if(l.nilg2 <= 0){
      std::stringstream ss, si;
      ss << l.nilg2;
      si << "NILG2[" << i << "]";
      printIO("LAMBDAS", si.str(), ss.str(), "> 0");
    }
    
    s.lambints.push_back(l);
    is >> dum;
    
  }
  e=dum;
  return is;
}

istringstream &operator>>(istringstream &is,ilocalelev &s){
  string e;
  s.found=1;
  is >> s.ntles >> s.ntlefr >> s.ntlefu >> s.nlegrd >> s.ntlesa
     >> s.cles >> s.wles >> s.rles >> e;
  if(s.ntles < 0 || s.ntles > 2){
    std::stringstream ss;
    ss << s.ntles;
    printIO("LOCALELEV", "NTLES", ss.str(), "0,1,2");
  }
  if(s.ntlefr < 0 || s.ntlefr > 1){
    std::stringstream ss;
    ss << s.ntlefr;
    printIO("LOCALELEV", "NTLEFR", ss.str(), "0,1");
  }
  if(s.ntlefu < 0 || s.ntlefu > 1){
    std::stringstream ss;
    ss << s.ntlefu;
    printIO("LOCALELEV", "NTLEFU", ss.str(), "0,1");
  }
  if(s.nlegrd <= 1){
    std::stringstream ss;
    ss << s.nlegrd;
    printIO("LOCALELEV", "NLEGRD", ss.str(), "> 1");
  }
  if(s.ntlesa < 0 || s.ntlesa > 2){
    std::stringstream ss;
    ss << s.ntlesa;
    printIO("LOCALELEV", "NTLESA", ss.str(), "0,1,2");
  }
  if(s.cles < 0.0){
    std::stringstream ss;
    ss << s.cles;
    printIO("LOCALELEV", "CLES", ss.str(), ">= 0");
  }
  if(s.wles < 0.0){
    std::stringstream ss;
    ss << s.wles;
    printIO("LOCALELEV", "WLES", ss.str(), ">= 0");
  }
  if(s.rles < 0.0){
    std::stringstream ss;
    ss << s.rles;
    printIO("LOCALELEV", "RLES", ss.str(), ">= 0");
  }

  return is;
}

istringstream &operator>>(istringstream &is,imultibath &s){
  string e, alg;
  double temp0, tau;
  int last, combath, irbath;
  s.found=1;
  is >> alg;
  std::transform(alg.begin(), alg.end(), alg.begin(), ::tolower);
  if(alg=="weak-coupling") s.algorithm=0;
  else if(alg=="nose-hoover") s.algorithm=1;
  else if(alg=="nose-hoover-chains"){
    s.algorithm=2;
    is >> s.num;
    if(s.num < 2){
      std::stringstream ss;
      ss << s.num;
      printIO("MULTIBATH", "NUM", ss.str(), "> 2");
    }
  }
  else{
    std::stringstream ss(alg);
    if( !(ss >> s.algorithm) || s.algorithm < 0 || s.algorithm > 2){
      printIO("MULTIBATH", "algorithm", alg, 
	      "weak-coupling(0), nose-hoover(1), nose-hoover-chains(2)");
    }
  }
  
  is >> s.nbaths;
  if(s.nbaths < 0){
    std::stringstream ss;
    ss << s.nbaths;
    printIO("MULTIBATH", "NBATHS", ss.str(), ">0");
  }
  for(int i=0;i<s.nbaths;++i){
    is >> temp0 >> tau;
    s.temp0.push_back(temp0);
    s.tau.push_back(tau);

    if(temp0 < 0.0){
      std::stringstream ss, si;
      si << "TEMP0[" << i+1 << "]";
      ss << temp0;
      printIO("MULTIBATH", si.str(), ss.str(), ">= 0.0");
    }
    if(tau < 0.0 && tau !=-1){
      std::stringstream ss, si;
      si << "TAU[" << i+1 << "]";
      ss << tau;
      printIO("MULTIBATH", si.str(), ss.str(), ">= 0.0, -1");
    }
  }
  
  is >> s.dofset;
  if(s.dofset < 0){
    std::stringstream ss;
    ss << s.dofset;
    printIO("MULTIBATH", "DOFSET", ss.str(), ">0");
  }
 
  for(int i=0;i<s.dofset;++i){
    is >> last >> combath >> irbath;
    s.last.push_back(last); 
    s.combath.push_back(combath); 
    s.irbath.push_back(irbath);
    if(last < 1){
      std::stringstream ss, si;
      si << "LAST[" << i+1 << "]";
      ss << last;
      printIO("MULTIBATH", si.str(), ss.str(), "> 1");
    }
    if(combath < 1){
      std::stringstream ss, si;
      si << "COMBATH[" << i+1 << "]";
      ss << combath;
      printIO("MULTIBATH", si.str(), ss.str(), "> 1");
    }
    if(irbath < 1){
      std::stringstream ss, si;
      si << "IRBATH[" << i+1 << "]";
      ss << irbath;
      printIO("MULTIBATH", si.str(), ss.str(), "> 1");
    }
    
  }
  is >> e;
  return is;
}
istringstream &operator>>(istringstream &is, imulticell &s){
  string e;
  s.found=1;
  is >> s.ntm >> s.ncella >> s.ncellb >> s.ncellc
     >> s.tolpx >> s.tolpv >> s.tolpf >> s.tolpfw >> e;
  if(s.ntm < 0 || s.ntm > 1){
    std::stringstream ss;
    ss << s.ntm;
    printIO("MULTICELL", "NTM", ss.str(), "0,1");
  }
  if(s.ncella < 1 ){
    std::stringstream ss;
    ss << s.ncella;
    printIO("MULTICELL", "NCELLA", ss.str(), ">= 1");
  }
  if(s.ncellb < 1 ){
    std::stringstream ss;
    ss << s.ncellb;
    printIO("MULTICELL", "NCELLB", ss.str(), ">= 1");
  }
  if(s.ncellc < 1 ){
    std::stringstream ss;
    ss << s.ncellc;
    printIO("MULTICELL", "NCELLC", ss.str(), ">= 1");
  }
  if(s.tolpx <= 0.0){
    std::stringstream ss;
    ss << s.tolpx;
    printIO("MULTICELL", "TOLPX", ss.str(), "> 0.0");
  }
  if(s.tolpv <= 0.0){
    std::stringstream ss;
    ss << s.tolpv;
    printIO("MULTICELL", "TOLPV", ss.str(), "> 0.0");
  }
  if(s.tolpf <= 0.0){
    std::stringstream ss;
    ss << s.tolpf;
    printIO("MULTICELL", "TOLPF", ss.str(), "> 0.0");
  }
  if(s.tolpfw <= 0.0){
    std::stringstream ss;
    ss << s.tolpfw;
    printIO("MULTICELL", "TOLPFW", ss.str(), "> 0.0");
  }
  
  return is;
}
istringstream &operator>>(istringstream &is,ineighbourlist &s){
  string e;
  s.found=1;
  is >> s.nmprpl >> s.nuprpl >> s.rcprpl >> s.grprpl >> s.nmtwpl;
  is >> s.nutwpl >> s.rstwpl >> s.rltwpl >> s.nuirin >> s.nusrin;
  is >> s.nmtwin >> s.rctwin >> s.ncgcen;
  if(s.nmprpl < 0 || s.nmprpl > 3){
    std::stringstream ss;
    ss << s.nmprpl;
    printIO("NEIGBOURLIST", "NMPRPL", ss.str(), "0,1,2,3");
  }
  if(s.nuprpl < 0 ){
    std::stringstream ss;
    ss << s.nuprpl;
    printIO("NEIGBOURLIST", "NMPRPL", ss.str(), ">= 0");
  }
  if(s.rcprpl < 0){
    std::stringstream ss;
    ss << s.rcprpl;
    printIO("NEIGBOURLIST", "RCPRPL", ss.str(), ">= 0.0");
  }
  if(s.grprpl < 0 || s.grprpl > 1.0){
    std::stringstream ss;
    ss << s.grprpl;
    printIO("NEIGBOURLIST", "GRPRPL", ss.str(), ">= 0.0 and <= 1.0");
  }
  if(s.nmtwpl < 0 || s.nmtwpl > 3){
    std::stringstream ss;
    ss << s.nmtwpl;
    printIO("NEIGBOURLIST", "NMTWPL", ss.str(), "0,1,2,3");
  }
  if(s.nutwpl < 0 ){
    std::stringstream ss;
    ss << s.nutwpl;
    printIO("NEIGBOURLIST", "NUTWPL", ss.str(), ">= 0");
  }
  if(s.rstwpl < 0){
    std::stringstream ss;
    ss << s.rstwpl;
    printIO("NEIGBOURLIST", "RSTWPL", ss.str(), ">= 0");
  }
  if(s.rltwpl < s.rstwpl){
    std::stringstream ss, si;
    ss << s.nmprpl;
    si << ">= RSTWPL (" << s.rstwpl << ")";
    printIO("NEIGBOURLIST", "RLTWPL", ss.str(), si.str());
  }
  if(s.nuirin < 0){
    std::stringstream ss;
    ss << s.nuirin;
    printIO("NEIGBOURLIST", "NUIRIN", ss.str(), ">= 0");
  }
  if(s.nusrin < 0){
    std::stringstream ss;
    ss << s.nusrin;
    printIO("NEIGBOURLIST", "NMPRPL", ss.str(), ">= 0");
  }
  if(s.nmtwin < 0 || s.nmtwin > 2){
    std::stringstream ss;
    ss << s.nmtwin;
    printIO("NEIGBOURLIST", "NMTWIN", ss.str(), "0,1,2");
  }
  if(s.rctwin < 0.0){
    std::stringstream ss;
    ss << s.rctwin;
    printIO("NEIGBOURLIST", "RCTWIN", ss.str(), ">= 0.0");
  }
  if(s.ncgcen < -2 ){
    std::stringstream ss;
    ss << s.nmprpl;
    printIO("NEIGBOURLIST", "NCGCEN", ss.str(), ">= -2");
  }
  
  return is;
}

istringstream &operator>>(istringstream &is,inonbonded &s){
  string e;
  s.found=1;
  is >> s.nlrele;
  is >> s.appak >> s.rcrf >> s.epsrf;
  is >> s.nshape >> s.ashape >> s.na2clc >> s.tola2 >> s.epsls;
  is >> s.nkx >> s.nky >> s.nkz >> s.kcut;
  is >> s.ngx >> s.ngy >> s.ngz >> s.nasord >> s.nfdord >> s.nalias 
     >> s.nspord;
  is >> s.nqeval >> s.faccur >> s.nrdgrd >> s.nwrgrd;
  is >> s.nlrlj >> s.slvdns >> e;

  if(s.nlrele < -4 || s.nlrele > 4){
    std::stringstream ss;
    ss << s.nlrele;
    printIO("NONBONDED", "NLRELE", ss.str(), "-4,-3,-2,-1,0,1,2,3,4");
  }
  if(s.appak < 0.0){
    std::stringstream ss;
    ss << s.appak;
    printIO("NONBONDED", "APPAK", ss.str(), ">= 0.0");
  }
  if(s.rcrf < 0.0){
    std::stringstream ss;
    ss << s.rcrf;
    printIO("NONBONDED", "RCRF", ss.str(), ">= 0.0");
  }
  if(s.epsrf < 1.0 && s.epsrf!=0.0){
    std::stringstream ss;
    ss << s.epsrf;
    printIO("NONBONDED", "EPSRF", ss.str(), "0.0 (off) or >= 1.0");
  }
  if(s.nshape < -1 || s.nshape > 10){
    std::stringstream ss;
    ss << s.nshape;
    printIO("NONBONDED", "NSHAPE", ss.str(), "-1 .. 10");
  }
  if(s.ashape <= 0.0){
    std::stringstream ss;
    ss << s.ashape;
    printIO("NONBONDED", "ASHAPE", ss.str(), "> 0.0");
  }
  if(s.na2clc < 0 || s.na2clc > 4){
    std::stringstream ss;
    ss << s.na2clc;
    printIO("NONBONDED", "NA2CLC", ss.str(), "0,1,2,3,4");
  }
  if(s.tola2 <= 0.0){
    std::stringstream ss;
    ss << s.tola2;
    printIO("NONBONDED", "TOLA2", ss.str(), "> 0.0");
  }
  if(s.epsls < 1.0 && s.epsls != 0.0){
    std::stringstream ss;
    ss << s.epsls;
    printIO("NONBONDED", "EPSLS", ss.str(), "0.0 (off) or >= 1.0");
  }
  if(s.nkx <= 0 ){
    std::stringstream ss;
    ss << s.nkx;
    printIO("NONBONDED", "NKX", ss.str(), "> 0");
  }
  if(s.nky <= 0){
    std::stringstream ss;
    ss << s.nky;
    printIO("NONBONDED", "NKY", ss.str(), "> 0");
  }
  if(s.nkz <= 0){
    std::stringstream ss;
    ss << s.nkz;
    printIO("NONBONDED", "NKZ", ss.str(), "> 0");
  }
  if(s.kcut <= 0.0){
    std::stringstream ss;
    ss << s.kcut;
    printIO("NONBONDED", "KCUT", ss.str(), "> 0.0");
  }
  if(s.ngx <= 0){
    std::stringstream ss;
    ss << s.ngx;
    printIO("NONBONDED", "NGX", ss.str(), "> 0");
  }
  if(s.ngy <= 0){
    std::stringstream ss;
    ss << s.ngy;
    printIO("NONBONDED", "NGY", ss.str(), "> 0");
  }
  if(s.ngz <= 0){
    std::stringstream ss;
    ss << s.ngz;
    printIO("NONBONDED", "NGZ", ss.str(), "> 0");
  }
  if(s.nasord < 1 || s.nasord > 5){
    std::stringstream ss;
    ss << s.nasord;
    printIO("NONBONDED", "NASORD", ss.str(), "1,2,3,4,5");
  }
  if(s.nfdord < 0 || s.nfdord > 5){
    std::stringstream ss;
    ss << s.nfdord;
    printIO("NONBONDED", "NFDORD", ss.str(), "0,1,2,3,4,5");
  }
  if(s.nalias <= 0){
    std::stringstream ss;
    ss << s.nalias;
    printIO("NONBONDED", "NALIAS", ss.str(), "> 0");
  }
  if(s.nqeval < 0 ){
    std::stringstream ss;
    ss << s.nqeval;
    printIO("NONBONDED", "NQEVAL", ss.str(), ">= 0");
  }
  if(s.faccur <= 0.0){
    std::stringstream ss;
    ss << s.faccur;
    printIO("NONBONDED", "FACCUR", ss.str(), "> 0");
  }
  if(s.nrdgrd < 0 || s.nrdgrd > 1){
    std::stringstream ss;
    ss << s.nrdgrd;
    printIO("NONBONDED", "NRDGRD", ss.str(), "0,1");
  }
  if(s.nwrgrd < 0 || s.nwrgrd > 1){
    std::stringstream ss;
    ss << s.nwrgrd;
    printIO("NONBONDED", "NWRGRD", ss.str(), "0,1");
  }
  if(s.nlrlj < 0 || s.nlrlj > 1){
    std::stringstream ss;
    ss << s.nlrlj;
    printIO("NONBONDED", "NLRLJ", ss.str(), "0,1");
  }
  if(s.slvdns <= 0.0){
    std::stringstream ss;
    ss << s.slvdns;
    printIO("NONBONDED", "SLVDNS", ss.str(), "> 0.0");
  }

  return is;
}

istringstream &operator>>(istringstream &is,ioveralltransrot &s){
  string e;
  s.found=1;
  is >> s.ncmtr >> s.ncmro >> s.cmamx >> s.cmamy >> s.cmamz >> e;
  if(s.ncmtr < 0 || s.ncmtr>1){
    std::stringstream ss;
    ss << s.ncmtr;
    printIO("OVERALLTRANSROT", "NCMTR", ss.str(), "0,1");
  }
  if(s.ncmro < 0 || s.ncmro>1){
    std::stringstream ss;
    ss << s.ncmro;
    printIO("OVERALLTRANSROT", "NCMRO", ss.str(), "0,1");
  }
  return is;
}

istringstream &operator>>(istringstream &is,ipairlist &s){
  string e, alg, siz, typ;
  s.found=1;
  is >> alg;
  std::transform(alg.begin(), alg.end(), alg.begin(), ::tolower);
  if(alg=="standard") s.algorithm=1;
  else if(alg=="grid") s.algorithm=2;
  else if(alg=="vgrid") s.algorithm=3;
  else{
    std::stringstream ss(alg);
    if(!(ss>>s.algorithm) || s.algorithm<1 || s.algorithm>3)
      printIO("PAIRLIST", "algorithm", alg, "standard(1), grid(2), vgrid(3)");
  }
  is >> s.nsnb >> s.rcutp >> s.rcutl >> siz;
  std::transform(siz.begin(), siz.end(), siz.begin(), ::tolower);
  s.size = siz;
  if(siz!="auto"){
    double dsize;
    std::stringstream ss(siz);
    if(!(ss>>dsize)|| dsize <=0.0)
      printIO("PAIRLIST", "SIZE", siz, "> 0.0, or auto");
  }
  if(s.nsnb <=0){
    std::stringstream ss;
    ss << s.nsnb;
    printIO("PAIRLIST", "NSNB", ss.str(), "> 0");
  }
  if(s.rcutl <=0.0){
    std::stringstream ss;
    ss << s.rcutl;
    printIO("PAIRLIST", "RCUTL", ss.str(), "> 0.0");
  }
  if(s.rcutp <=0.0){
    std::stringstream ss;
    ss << s.rcutp;
    printIO("PAIRLIST", "RCUTP", ss.str(), "> 0.0");
  }

  is >> typ;
  std::transform(typ.begin(), typ.end(), typ.begin(), ::tolower);
  if(typ=="chargegroup") s.type=0;
  else if(typ=="atomic") s.type=1;
  else {
    std::stringstream ss(typ);
    if(!(ss>>s.type) || s.type < 0 || s.type > 1)
      printIO("PAIRLIST", "TYPE", typ, "chargegroup(0), atomic(1)");
  }
  is >> e;
  return is;
}

istringstream &operator>>(istringstream &is,ipathint &s){
  string e;
  s.found=1;
  is >> s.ntpi >> e;
  if(s.ntpi<0 || s.ntpi > 1){
    std::stringstream ss;
    ss << s.ntpi;
    printIO("PATHINT", "NTPI", ss.str(), "0,1");
  }
  return is;
}

istringstream &operator>>(istringstream &is,iperscale &s){
  string e, dum;
  s.found=1;
  is >> dum;
  std::transform(dum.begin(), dum.end(), dum.begin(), ::tolower); 
  s.restype=dum;
  if(dum!="jrest")
    printIO("PERSCALE", "RESTYPE", dum, "jrest");

  is >> s.kdih >> s.kj >> s.t >> s.diff >> s.ratio >> s.read >> e;
  if(s.kdih < 0.0){
    std::stringstream ss;
    ss << s.kdih;
    printIO("PERSCALE", "KDIH", ss.str(), ">= 0.0");
  }
  if(s.kj < 0.0){
    std::stringstream ss;
    ss << s.kj;
    printIO("PERSCALE", "KJ", ss.str(), ">= 0.0");
  }
  if(s.t <= 0){
    std::stringstream ss;
    ss << s.t;
    printIO("PERSCALE", "T", ss.str(), "> 0");
  }
  if(s.diff < 0.0){
    std::stringstream ss;
    ss << s.diff;
    printIO("PERSCALE", "DIFF", ss.str(), ">= 0.0");
  }
  if(s.ratio <= 0.0){
    std::stringstream ss;
    ss << s.ratio;
    printIO("PERSCALE", "RATIO", ss.str(), "> 0.0");
  }
  if(s.read < 0 || s.read > 1){
    std::stringstream ss;
    ss << s.read;
    printIO("PERSCALE", "READ", ss.str(), "0,1");
  }
  return is;
}

istringstream &operator>>(istringstream &is,iperturbation &s){
  string e;
  s.found=1;
  is >> s.ntg >> s.nrdgl >> s.rlam >> s.dlamt
     >> s.alphlj >> s.alphc >> s.nlam >> s.nscale >> e;
  if(s.ntg < 0 || s.ntg > 1){
    std::stringstream ss;
    ss << s.ntg;
    printIO("PERTURBATION", "NTG", ss.str(), "0,1");
  }
  if(s.nrdgl < 0 || s.nrdgl > 1){
    std::stringstream ss;
    ss << s.ntg;
    printIO("PERTURBATION", "NRDGL", ss.str(), "0,1");
  }
  if(s.rlam < 0.0 || s.rlam > 1.0){
    std::stringstream ss;
    ss << s.rlam;
    printIO("PERTURBATION", "RLAM", ss.str(), ">= 0.0 and <= 1.0");
  }
  if(s.dlamt < 0.0){
    std::stringstream ss;
    ss << s.dlamt;
    printIO("PERTURBATION", "DLAMT", ss.str(), ">= 0.0");
  }
  if(s.alphlj < 0.0){
    std::stringstream ss;
    ss << s.alphlj;
    printIO("PERTURBATION", "ALPHLJ", ss.str(), ">= 0.0");
  }
  if(s.alphc < 0.0){
    std::stringstream ss;
    ss << s.alphc;
    printIO("PERTURBATION", "ALPHC", ss.str(), ">= 0.0");
  }
  if(s.nlam < 0){
    std::stringstream ss;
    ss << s.nlam;
    printIO("PERTURBATION", "NLAM", ss.str(), "> 0");
  }
  if(s.nscale < 0 || s.nscale > 2){
    std::stringstream ss;
    ss << s.nscale;
    printIO("PERTURBATION", "NSCALE", ss.str(), "0,1,2");
  }
  return is;
}

istringstream &operator>>(istringstream &is, ipolarise &s){
  string e;
  s.found=1;
  is >> s.cos >> s.efield >> s.minfield >> s.damp >> s.write >> e;
  if(s.cos < 0 || s.cos > 1){
    std::stringstream ss;
    ss << s.cos;
    printIO("POLARISE", "COS", ss.str(), "0,1");
  }
  if(s.efield < 0 || s.efield > 1){
    std::stringstream ss;
    ss << s.efield;
    printIO("POLARISE", "EFIELD", ss.str(), "0,1");
  }
  if(s.minfield <= 0.0){
    std::stringstream ss;
    ss << s.minfield;
    printIO("POLARISE", "MINFIELD", ss.str(), "> 0.0");
  }
  if(s.damp < 0 || s.damp > 1){
    std::stringstream ss;
    ss << s.damp;
    printIO("POLARISE", "DAMP", ss.str(), "0,1");
  }
  if(s.write < 0){
    std::stringstream ss;
    ss << s.write;
    printIO("POLARISE", "WRITE", ss.str(), "> 0");
  }
     
  return is;
}

istringstream &operator>>(istringstream &is,ipositionres &s){
  string e;
  s.found=1;
  is >> s.ntpor >> s.ntporb >> s.ntpors >> s.cpor >> e;
  if(s.ntpor < 0 || s.ntpor > 3){
    std::stringstream ss;
    ss << s.ntpor;
    printIO("POSITIONRES", "NTPOR", ss.str(), "0,1,2,3");
  }
  if(s.ntporb < 0 || s.ntporb > 1){
    std::stringstream ss;
    ss << s.ntporb;
    printIO("POSITIONRES", "NTPORB", ss.str(), "0,1");
  }
  if(s.ntpors < 0 || s.ntpors > 1){
    std::stringstream ss;
    ss << s.ntpors;
    printIO("POSITIONRES", "NTPORS", ss.str(), "0,1");
  }
  if(s.cpor < 0){
    std::stringstream ss;
    ss << s.cpor;
    printIO("POSITIONRES", "CPOR", ss.str(), ">= 0.0");
  }
  
  return is;
}

istringstream &operator>>(istringstream &is,ipressurescale &s){
  string e, coup, sca, vir;
  s.found=1;

  is >> coup;
  std::transform(coup.begin(), coup.end(), coup.begin(), ::tolower);
  if(coup=="off") s.couple=0;
  else if(coup=="calc") s.couple=1;
  else if(coup=="scale") s.couple=2;
  else{
    std::stringstream ss(coup);
    if(!(ss >> s.couple) || s.couple < 0 || s.couple > 2)
      printIO("PRESSURESCALE", "COUPLE", coup, "off(0),calc(1),scale(2)");
  }
  is >> sca;
  std::transform(sca.begin(), sca.end(), sca.begin(), ::tolower);
  if(sca=="off") s.scale=0;
  else if(sca=="iso") s.scale=1;
  else if(sca=="aniso") s.scale=2;
  else if(sca=="full") s.scale=3;
  else{
    std::stringstream ss(sca);
    if(!(ss >> s.scale) || s.scale < 0 || s.scale > 3)
      printIO("PRESSURESCALE", "SCALE", sca, "off(0),iso(1),aniso(2),full(3)");
  }
  is >>  s.comp >> s.taup >> vir;
  std::transform(vir.begin(), vir.end(), vir.begin(), ::tolower);
  if(s.comp <= 0.0){
    std::stringstream ss;
    ss << s.comp;
    printIO("PRESSURESCALE", "COMP", ss.str(), "> 0.0");
  }
  if(s.taup <=0.0){
    std::stringstream ss;
    ss << s.taup;
    printIO("PRESSURESCALE", "COMP", ss.str(), "> 0.0");
  }
  if(vir=="none") s.virial=0;
  else if(vir=="atomic") s.virial=1;
  else if(vir=="molecular") s.virial=2;
  else{
    std::stringstream ss(vir);
    if(!(ss >> s.virial) || s.virial < 0 || s.virial > 2)
      printIO("PRESSURESCALE", "VIRIAL", vir, "none(0),atomic(1),molecular(2)");
  }
  
  is >> s.pres0[0][0] >> s.pres0[0][1] >> s.pres0[0][2]
     >> s.pres0[1][0] >> s.pres0[1][1] >> s.pres0[1][2]
     >> s.pres0[2][0] >> s.pres0[2][1] >> s.pres0[2][2] >> e;

  return is;
}

istringstream &operator>>(istringstream &is,iprintout &s){
  string e;
  s.found=1;
  is >> s.ntpr >> s.ntpp >> e;
  if(s.ntpr < 0){
    std::stringstream ss;
    ss << s.ntpr;
    printIO("PRINTOUT", "NTPR", ss.str(), ">= 0");
  }
  if(s.ntpp < 0 || s.ntpp > 1){
    std::stringstream ss;
    ss << s.ntpp;
    printIO("PRINTOUT", "NTPP", ss.str(), "0,1");
  }
  return is;
}

istringstream &operator>>(istringstream &is, irandomnumbers &s){
  string e;
  s.found=1;
  is >> s.ntrng >> s.ntgsl >> e;
  if(s.ntrng < 0 || s.ntrng > 1){
    std::stringstream ss;
    ss << s.ntrng;
    printIO("RANDOMNUMBERS", "NTRNG", ss.str(), "0,1");
  }
  if(s.ntgsl < -1){
    std::stringstream ss;
    ss << s.ntgsl;
    printIO("RANDOMNUMBERS", "NTGSL", ss.str(), "> -1");
  }
  return is;
}

istringstream &operator>>(istringstream &is, ireadtraj &s){
  string e;
  s.found=1;
  is >> s.ntrd >> s.ntrn >> s.ntrb >> s.ntshk >> e;
  if(s.ntrd < 0 || s.ntrd > 1){
    std::stringstream ss;
    ss << s.ntrd;
    printIO("READTRAJ", "NTRD", ss.str(), "0,1");
  }
  if(s.ntrn < 1 || s.ntrn > 18){
    std::stringstream ss;
    ss << s.ntrn;
    printIO("READTRAJ", "NTRN", ss.str(), "1..18");
  }
  if(s.ntrb < 0 || s.ntrd > 1){
    std::stringstream ss;
    ss << s.ntrb;
    printIO("READTRAJ", "NTRB", ss.str(), "0,1");
  }
  if(s.ntshk < 0 || s.ntshk > 1){
    std::stringstream ss;
    ss << s.ntshk;
    printIO("READTRAJ", "NTSHK", ss.str(), "0,1");
  }
  
  return is;
}

istringstream &operator>>(istringstream &is,ireplica &s){
  string e;
  int nret, nrelam;
  double dum;
  s.found=1;
  is >> nret;
  if(nret< 0){
    std::stringstream ss;
    ss << nret;
    printIO("REPLICA", "NRET", ss.str(), ">= 0");
  }
  for(int i=0; i<nret; ++i){
    is >> dum;
    s.ret.push_back(dum);
    if(dum < 0){
      std::stringstream ss, si;
      ss << dum;
      si << "RET[" << i+1 << "]";
      printIO("REPLICA", si.str(), ss.str(), ">= 0.0");
    }
  }
  is >> s.lrescale;
  if(s.lrescale < 0 || s.lrescale > 1){
    std::stringstream ss;
    ss << s.lrescale;
    printIO("REPLICA", "LRESCALE", ss.str(), "0,1");
  }
  is >> nrelam;
  if(nrelam< 0){
    std::stringstream ss;
    ss << nrelam;
    printIO("REPLICA", "NRELAM", ss.str(), ">= 0");
  }
  for(int i=0; i<nrelam; ++i){
    is >> dum;
    s.relam.push_back(dum);
    if(dum < 0){
      std::stringstream ss, si;
      ss << dum;
      si << "RELAM[" << i+1 << "]";
      printIO("REPLICA", si.str(), ss.str(), ">= 0.0");
    }
    is >> dum;
    s.rets.push_back(dum);
    if(dum < 0){
      std::stringstream ss, si;
      ss << dum;
      si << "RETS[" << i+1 << "]";
      printIO("REPLICA", si.str(), ss.str(), ">= 0.0");
    }
  }
  is >> s.nretrial >> s.nrequil >> s.nrejob >> s.nrewrt >> e;
  if(s.nretrial< 0){
    std::stringstream ss;
    ss << s.nretrial;
    printIO("REPLICA", "NRETRIAL", ss.str(), ">= 0");
  }
  if(s.nrequil< 0){
    std::stringstream ss;
    ss << s.nrequil;
    printIO("REPLICA", "NREQUIL", ss.str(), ">= 0");
  }
  if(s.nrejob <= 0){
    std::stringstream ss;
    ss << s.nrejob;
    printIO("REPLICA", "NREJOB", ss.str(), "> 0");
  }
  if(s.nrewrt< 0){
    std::stringstream ss;
    ss << s.nrewrt;
    printIO("REPLICA", "NREWRT", ss.str(), ">= 0");
  }
  
  return is;
}
istringstream &operator>>(istringstream &is,irottrans &s){
  string e;
  s.found=1;
  is >> s.rtc >> s.rtclast >> e;
  if(s.rtc < 0 || s.rtc > 1){
    std::stringstream ss;
    ss << s.rtc;
    printIO("ROTTRANS", "RTC", ss.str(), "0,1");
  }
  if(s.rtclast <= 0){
    std::stringstream ss;
    ss << s.rtclast;
    printIO("ROTTRANS", "RTCLAST", ss.str(), "> 0");
  }
  
  return is;
}

istringstream &operator>>(istringstream &is, istep &s){
  string e;
  s.found=1;
  is >> s.nstlim >> s.t >> s.dt >> e;
  if(s.nstlim < 0){
    std::stringstream ss;
    ss << s.nstlim;
    printIO("STEP", "NSTLIM", ss.str(), ">= 0");
  }
  if(s.t < 0.0){
    std::stringstream ss;
    ss << s.t;
    printIO("STEP", "T", ss.str(), ">= 0.0");
  }
  if(s.dt <= 0){
    std::stringstream ss;
    ss << s.dt;
    printIO("STEP", "DT", ss.str(), "> 0.0");
  }

  return is;
}

istringstream &operator>>(istringstream &is, istochdyn &s){
  string e;
  s.found=1;
  is >> s.ntsd >> s.ntfr >> s.nsfr >> s.nbref >> s.rcutf >> s.cfric >> s.tempsd >> e;
  if(s.ntsd < 0 || s.ntsd > 1){
    std::stringstream ss;
    ss << s.ntsd;
    printIO("STOCHDYN", "NTSD", ss.str(), "0,1");
  }
  if(s.ntfr < 0 || s.ntfr > 3){
    std::stringstream ss;
    ss << s.ntfr;
    printIO("STOCHDYN", "NTFR", ss.str(), "0,1,2,3");
  }
  if(s.nsfr <= 0){
    std::stringstream ss;
    ss << s.nsfr;
    printIO("STOCHDYN", "NSFR", ss.str(), "> 0");
  }
  if(s.nbref <= 0){
    std::stringstream ss;
    ss << s.nbref;
    printIO("STOCHDYN", "NBREF", ss.str(), "> 0");
  }
  if(s.rcutf < 0 ){
    std::stringstream ss;
    ss << s.rcutf;
    printIO("STOCHDYN", "RCUTF", ss.str(), ">= 0.0");
  }
  if(s.cfric < 0){
    std::stringstream ss;
    ss << s.cfric;
    printIO("STOCHDYN", "CFRIC", ss.str(), ">= 0.0");
  }
  if(s.tempsd < 0){
    std::stringstream ss;
    ss << s.tempsd;
    printIO("STOCHDYN", "TEMPSD", ss.str(), ">= 0.0");
  }

  return is;
}

istringstream &operator>>(istringstream &is, isystem &s){
  string e;
  s.found=1;
  is >> s.npm >> s.nsm >> e;
  if(s.npm < 0){
    std::stringstream ss;
    ss << s.npm;
    printIO("SYSTEM", "NPM", ss.str(), ">= 0");
  }
  if(s.nsm < 0){
    std::stringstream ss;
    ss << s.nsm;
    printIO("SYSTEM", "NSM", ss.str(), ">= 0");
  }

  return is;
}

istringstream &operator>>(istringstream &is,ithermostat &s){
  string e;
  int idum, dum, ntbth, ntset;
  double tau;
  s.found=1;
  is >> s.ntt >> ntbth >> ntset;
  if(s.ntt < 0 || s.ntt >1){
    std::stringstream ss;
    ss << s.ntt;
    printIO("THERMOSTAT", "NTT", ss.str(), "0,1");
  }
  if(ntbth < 0){
    std::stringstream ss;
    ss << ntbth;
    printIO("THERMOSTAT", "NTBTH", ss.str(), ">=0");
  }
  if(ntset < 0 ){
    std::stringstream ss;
    ss << ntset;
    printIO("THERMOSTAT", "NTSET", ss.str(), ">= 0");
  }
  for(int i=0;i<ntbth;++i){
    class ithermostat::tbath b;
    is >> idum >> b.ntbtyp>> b.tembth >> dum;
    if(idum != i+1){
      std::stringstream ss, si, sj;
      ss << idum;
      si << "I[" << i+1 << "]";
      sj << i+1;
      printIO("THERMOSTAT", si.str(), ss.str(), sj.str());
    }
    if(b.ntbtyp < 0 || b.ntbtyp > 3){
      std::stringstream ss, si;
      ss << b.ntbtyp;
      si << "NTBTYP[" << i+1 << "]";
      printIO("THERMOSTAT", si.str(), ss.str(), "0,1,2,3");
    }
    if(b.tembth < 0.0){
      std::stringstream ss, si;
      ss << b.tembth;
      si << "TEMBTH[" << i+1 << "]";
      printIO("THERMOSTAT", si.str(), ss.str(), ">= 0.0");
    }
    if(dum < 0){
      std::stringstream ss, si;
      ss << dum;
      si << "NTBVAR[" << i+1 << "]";
      printIO("THERMOSTAT", si.str(), ss.str(), ">= 0.0");
    }
    
    for(int j=0;j<dum; ++j){
       is >> tau;
       b.taubth.push_back(tau);
       if(tau < 0.0){
	 std::stringstream ss, si;
	 ss << tau;
	 si << "TAUBTH[" << i+1 << "," << j+1 << "]";
	 printIO("THERMOSTAT", si.str(), ss.str(), ">= 0.0");
       }
    }
    s.baths.push_back(b);
  }
  
  for(int i=0;i<ntset;++i){
    class ithermostat::dofgroup d;
    is >> d.ntscpl >> d.ntstyp >> d.ntscns >> d.ntsgt;
    if(d.ntscpl < 0 || d.ntscpl > ntbth){
      std::stringstream ss, si, sj;
      ss << d.ntscpl;
      si << "NTSCPL[" << i+1 << "]";
      sj << "0..NTBTH (" << ntbth << ")";
      printIO("THERMOSTAT", si.str(), ss.str(), sj.str());
    }
    if(d.ntstyp < 0 || d.ntstyp > 2){
      std::stringstream ss, si;
      ss << d.ntstyp;
      si << "NTSTYP[" << i+1 << "]";
      printIO("THERMOSTAT", si.str(), ss.str(), "0,1,2");
    }
    if(d.ntscns < 0 || d.ntscns > 1){
      std::stringstream ss, si;
      ss << d.ntscns;
      si << "NTSCNS[" << i+1 << "]";
      printIO("THERMOSTAT", si.str(), ss.str(), "0,1");
    }
    if(d.ntsgt < -2){
      std::stringstream ss, si;
      ss << d.ntsgt;
      si << "NTSGT[" << i+1 << "]";
      printIO("THERMOSTAT", si.str(), ss.str(), ">= -2");
    }
    int gt=0, ntsgtg;
    if (d.ntsgt <=0) gt=0;
    else gt=d.ntsgt;
    for(int k=0; k< gt; k++){
      is >> ntsgtg;
      d.ntsgtg.push_back(ntsgtg);
      if(ntsgtg <= 0){
	std::stringstream ss, si;
	ss << ntsgtg;
	si << "NTSGTG[" << i+1 << "," << k+1 << "]";
	printIO("THERMOSTAT", si.str(), ss.str(), ">= 1");
      }
    }
    s.dofgroups.push_back(d);
  }
  
  is >> e;
  return is;
}

istringstream &operator>>(istringstream &is,iumbrella &s){
  string e;
  s.found=1;
  is >> s.ntus >> s.uscst1 >> s.uscst2 >> s.usref1 >> s.usref2 >> e;
  if(s.ntus < 0 || s.ntus > 1){
    std::stringstream ss;
    ss << s.ntus;
    printIO("UMBRELLA", "NTUS", ss.str(), "0,1");
  }
  if(s.uscst1 < 0.0){
    std::stringstream ss;
    ss << s.uscst1;
    printIO("UMBRELLA", "USCST1", ss.str(), ">= 0.0");
  }
  if(s.uscst2 < 0.0){
    std::stringstream ss;
    ss << s.uscst2;
    printIO("UMBRELLA", "USCST2", ss.str(), ">= 0.0");
  }
  if(s.usref1 < 0.0){
    std::stringstream ss;
    ss << s.usref1;
    printIO("UMBRELLA", "USREF1", ss.str(), ">= 0.0");
  }
  if(s.usref2 < 0.0){
    std::stringstream ss;
    ss << s.usref2;
    printIO("UMBRELLA", "USREF2", ss.str(), ">= 0.0");
  }
  
  return is;
}

istringstream &operator>>(istringstream &is,ivirial &s){
  string e;
  s.found=1;
  is >> s.ntv >> s.ntvg >> e;
  if(s.ntv < 0 || s.ntv > 1){
    std::stringstream ss;
    ss << s.ntv;
    printIO("VIRIAL", "NTV", ss.str(), "0,1");
  }
  if(s.ntvg < 0 || s.ntvg > 1){
    std::stringstream ss;
    ss << s.ntvg;
    printIO("VIRIAL", "NTVG", ss.str(), "0,1,2,3");
  }

  return is;
}

istringstream &operator>>(istringstream &is,iwritetraj &s){
  string e;
  s.found=1;
  is >> s.ntwx >> s.ntwse >> s.ntwv >> s.ntwf >> s.ntwe >> s.ntwg >> s.ntwb >> e;
  if(s.ntwse < 0.0){
    std::stringstream ss;
    ss << s.ntwse;
    printIO("WRITETRAJ", "NTWSE", ss.str(), ">= 0.0");
  }
  if(s.ntwf < 0.0){
    std::stringstream ss;
    ss << s.ntwf;
    printIO("WRITETRAJ", "NTWF", ss.str(), ">= 0.0");
  }
  if(s.ntwe < 0.0){
    std::stringstream ss;
    ss << s.ntwe;
    printIO("WRITETRAJ", "NTWE", ss.str(), ">= 0.0");
  }
  if(s.ntwg < 0.0){
    std::stringstream ss;
    ss << s.ntwg;
    printIO("WRITETRAJ", "NTWG", ss.str(), ">= 0.0");
  }
  if(s.ntwb < 0.0){
    std::stringstream ss;
    ss << s.ntwb;
    printIO("WRITETRAJ", "NTWB", ss.str(), ">= 0.0");
  }
  
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
	case barostatblock:         bfstream >> gin.barostat;         break;
	case boundcondblock:        bfstream >> gin.boundcond;        break;
	case cgrainblock:           bfstream >> gin.cgrain;           break;
	case comtransrotblock:      bfstream >> gin.comtransrot;      break;
	case consistencycheckblock: bfstream >> gin.consistencycheck; break;
	case constraintblock:       bfstream >> gin.constraint;       break;
	case covalentformblock:     bfstream >> gin.covalentform;     break;
	case debugblock:            bfstream >> gin.debug;            break;
	case dihedralresblock:      bfstream >> gin.dihedralres;      break;
	case distanceresblock:      bfstream >> gin.distanceres;      break;
	case energyminblock:        bfstream >> gin.energymin;        break;
	case ewarnblock:            bfstream >> gin.ewarn;            break;
	case forceblock:            bfstream >> gin.force;            break;
	case geomconstraintsblock:  bfstream >> gin.geomconstraints;  break;
	case gromos96compatblock:   bfstream >> gin.gromos96compat;   break;
	case initialiseblock:       bfstream >> gin.initialise;       break;
	case innerloopblock:        bfstream >> gin.innerloop;        break;
	case integrateblock:        bfstream >> gin.integrate;        break;
	case jvalueresblock:        bfstream >> gin.jvalueres;        break;
	case lambdasblock:          bfstream >> gin.lambdas;          break;
	case localelevblock:        bfstream >> gin.localelev;        break;
	case multibathblock:        bfstream >> gin.multibath;        break;
	case multicellblock:        bfstream >> gin.multicell;        break;
	case neighbourlistblock:    bfstream >> gin.neighbourlist;    break;
	case nonbondedblock:        bfstream >> gin.nonbonded;        break;
	case overalltransrotblock:  bfstream >> gin.overalltransrot;  break;
	case pairlistblock:         bfstream >> gin.pairlist;         break;
	case pathintblock:          bfstream >> gin.pathint;          break;
	case perscaleblock:         bfstream >> gin.perscale;         break;
	case perturbationblock:     bfstream >> gin.perturbation;     break;
	case polariseblock:         bfstream >> gin.polarise;         break;
	case positionresblock:      bfstream >> gin.positionres;      break;
	case pressurescaleblock:    bfstream >> gin.pressurescale;    break;
	case printoutblock:         bfstream >> gin.printout;         break;
	case randomnumbersblock:    bfstream >> gin.randomnumbers;    break;
	case readtrajblock:         bfstream >> gin.readtraj;         break;
	case replicablock:          bfstream >> gin.replica;          break;
	case rottransblock:         bfstream >> gin.rottrans;         break;
	case stepblock:             bfstream >> gin.step;             break;
	case stochdynblock:         bfstream >> gin.stochdyn;         break;
	case systemblock:           bfstream >> gin.system;           break;
	case thermostatblock:       bfstream >> gin.thermostat;       break;
	case umbrellablock:         bfstream >> gin.umbrella;         break;
	case virialblock:           bfstream >> gin.virial;           break;
	case writetrajblock:        bfstream >> gin.writetraj;        break;
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
      // use the gio read_box functions?
      istringstream iss(buffer[0]);
      iss >> s.box[0] >> s.box[1] >> s.box[2] >> e;
    }
    if (first=="TRICLINICBOX"){

      int ntb;
      gmath::Vec k,l,m;
      istringstream iss(buffer[0]);

      iss >> ntb;
      s.box.boxformat()=gcore::Box::triclinicbox;
      
      s.box.setNtb(gcore::Box::boxshape_enum(ntb));
      
      iss.str(buffer[1]);
      iss >> k[0] >> l[0] >> m[0];
      iss.str(buffer[2]);
      iss >> k[1] >> l[1] >> m[1];
      iss.str(buffer[3]);
      iss >> k[2] >> l[2] >> m[2];
      s.box.K()=k;
      s.box.L()=l;
      s.box.M()=m;
    }
    if(first=="GENBOX"){
      s.box.boxformat()=gcore::Box::genbox;
      int ntb;
      gmath::Vec k,l,m;
      istringstream iss(buffer[0]);
      iss >> ntb;
      s.box.setNtb(gcore::Box::boxshape_enum(ntb));
      // I don't trust the gio reading of GENBOX.
      // so here, I just store the genbox data in the klm vectors.
      // that is very ugly. But we know what the format is. 
      // We only need it here to calculate the boxsize
      iss.str(buffer[1]);
      iss >> k[0] >> l[0] >> m[0];
      iss.str(buffer[2]);
      iss >> k[1] >> l[1] >> m[1];
      iss.str(buffer[3]);
      iss >> k[2] >> l[2] >> m[2];
      s.box.K()=k;
      s.box.L()=l;
      s.box.M()=m;

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
	case start_timetemplate: os << d_time+number*d_dt;     break;
	case end_timetemplate:   os << d_time+(number+1)*d_dt; break;
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
  // MOLECULAR SYSTEM

  // SYSTEM (g96, promd, md++)
  if(gin.system.found)
    os << "SYSTEM\n"
       << "#      NPM      NSM\n"
       << setw(10) << gin.system.npm
       << setw(9)  << gin.system.nsm
       << "\nEND\n";

  // METHOD EMPLOYED

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
  // STEP (promd, md++, g96)
  if(gin.step.found){
    os << "STEP\n"
       << "#   NSTLIM         T        DT\n"
       << setw(10) << gin.step.nstlim
       << setw(10) << gin.step.t
       << setw(10) << gin.step.dt
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

  // SPACIAL BOUNDARY CONDISTIONS

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

  // THERMODYNAMIC BOUNDARY CONDITIONS
  // THERMOSTAT (promd)
  if(gin.thermostat.found){
    os << "THERMOSTAT\n"
       << "#      NTT     NTBTH     NTSET\n"
       << setw(10) << gin.thermostat.ntt
       << setw(10) << gin.thermostat.baths.size()
       << setw(10) << gin.thermostat.dofgroups.size()
       << "\n"
       << "# I = 1 ... NTBTH\n"
       << "# I TEMBTH(I)  NTBVAR(I)   TAUBTH(I,1...NTVAR)\n";
    for(unsigned int i=0; i<gin.thermostat.baths.size(); ++i){
      os << setw(3)  << i+1
         << setw(10) << gin.thermostat.baths[i].tembth
	 << setw(10) << gin.thermostat.baths[i].taubth.size();
      
      for(unsigned int j=0; j<gin.thermostat.baths[i].taubth.size(); ++j){
        os << setw(7) << gin.thermostat.baths[i].taubth[j];
      }
      os << "\n";
    }
    os << "# J = 1 ... NTGRP\n"
       << "# NTSCPL(J) NTSTYP(J) NTSCNS(J)  NTSGT(J)  NTSGTG(J,1...NTGT(J))\n";
    for(unsigned int j=0; j<gin.thermostat.dofgroups.size(); ++j){
      os << setw(10) << gin.thermostat.dofgroups[j].ntscpl
         << setw(10) << gin.thermostat.dofgroups[j].ntstyp
         << setw(10) << gin.thermostat.dofgroups[j].ntscns
	 << setw(10) << gin.thermostat.dofgroups[j].ntsgt;
      if(gin.thermostat.dofgroups[j].ntsgt>0){
        for(unsigned int k=0;
	    k<gin.thermostat.dofgroups[j].ntsgtg.size(); ++k){
          os << setw(7) << gin.thermostat.dofgroups[j].ntsgtg[k];
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
    if(gin.multibath.algorithm==2){
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
    os << "# TEMP0(1 ... NBATHS)  TAU(1 ... NBATHS)\n";
    for(int i=0; i<gin.multibath.nbaths; ++i){
      os << setw(10) << gin.multibath.temp0[i]
	 << setw(10) << gin.multibath.tau[i];
    }
    os << "\n";
    os << "#   DOFSET: number of distiguishable sets of d.o.f.\n";
    os << setw(10) << gin.multibath.dofset << "\n";
    os << "# LAST(1 ... DOFSET)  COMBATH(1 ... DOFSET)  IRBATH(1 ... DOFSET)\n";
    for(int i=0; i<gin.multibath.dofset; ++i){
      os << setw(10) << gin.multibath.last[i]
	 << setw(10) << gin.multibath.combath[i]
	 << setw(10) << gin.multibath.irbath[i];
    }
    os << "\nEND\n";
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

  // INTERACTION EVALUATION

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
  // GEOMCONSTRAINTS (promd)
  if(gin.geomconstraints.found){
    os << "GEOMCONSTRAINTS\n"
       << "#    NTCPH     NTCPN      NTCS    SHKTOL\n"
       << setw(10) << gin.geomconstraints.ntcph
       << setw(10) << gin.geomconstraints.ntcpn
       << setw(10) << gin.geomconstraints.ntcs
       << setw(10) << gin.geomconstraints.shktol
       << "\nEND\n";
  }
  // CONSTRAINT (md++)
  if(gin.constraint.found){
    os << "CONSTRAINT\n"
       << "# NTC\n"
       << setw(5) << gin.constraint.ntc
       << "\n";
    if(gin.constraint.ntcp == 3){
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
    if(gin.constraint.ntcs != 3){
      os << "#      NTCS   NTCS0(1)\n"
         << setw(11) << gin.constraint.ntcs
         << setw(10) << gin.constraint.ntcs0[0];
    }
    else{
      os << "#      NTCS\n"
         << setw(11) << gin.constraint.ntcs;
    }
    os << "\nEND\n";
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
  // POLARISE 
  if(gin.polarise.found){
    os << "POLARISE\n"
       << "#     COS    EFIELD      DAMP     WRITE\n"
       << setw(10) << gin.polarise.cos
       << setw(10) << gin.polarise.efield
       << setw(10) << gin.polarise.damp
       << setw(10) << gin.polarise.write
       << "\nEND\n";
  }
  // INTEGRATE (md++)
  if(gin.integrate.found){
    os << "INTEGRATE\n"
       << "#   NINT\n"
       << setw(8) << gin.integrate.nint
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
  // ROTTRANS (md++)
  if(gin.rottrans.found){
    os << "ROTTRANS\n"
       << "#      RTC   RTCLAST\n"
       << setw(10) << gin.rottrans.rtc
       << setw(10) << gin.rottrans.rtclast
       << "\nEND\n";
  }
  // INNERLOOP (md++)
  if(gin.innerloop.found){
    os << "INNERLOOP\n"
       << "#     NTILM      NTILS\n"
       << setw(10) << gin.innerloop.ntilm
       << setw(10) << gin.innerloop.ntils
       << "\nEND\n";
  }

  // PAIRLIST GENERATION

  // NEIGHBOURLIST (promd)
  if(gin.neighbourlist.found){
    os << "NEIGHBOURLIST\n"
       << "# NMPRPL   NUPRPL   RCPRPL   GRPRPL\n"
       << setw(7) << gin.neighbourlist.nmprpl
       << setw(7) << gin.neighbourlist.nuprpl
       << setw(7) << gin.neighbourlist.rcprpl
       << setw(7) << gin.neighbourlist.grprpl
       << "\n"
       << "# NMTWPL   NUTWPL   RSTWPL:   TLTWPL\n"
       << setw(7) << gin.neighbourlist.nmtwpl
       << setw(7) << gin.neighbourlist.nutwpl
       << setw(7) << gin.neighbourlist.rstwpl
       << setw(7) << gin.neighbourlist.rltwpl
       << "\n"
       << "# NUIRIN   NUSRIN\n"
       << setw(7) << gin.neighbourlist.nuirin
       << setw(7) << gin.neighbourlist.nusrin
       << "\n"
       << "# NMTWIN   RCTWIN\n"
       << setw(7) << gin.neighbourlist.nmtwin
       << setw(7) << gin.neighbourlist.rctwin
       << "\n"
       << "# NCGCEN\n"
       << setw(7) << gin.neighbourlist.ncgcen
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

  // LONGRANGE INTERACTIONS

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

  // INITIALISATION OF THE RUN

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
  // RANDOMNUMBERS (md++)
  if(gin.randomnumbers.found){
    os << "RANDOMNUMBERS\n"
       << "#   NTRNG   NTGSL\n"
       << setw(8) << gin.randomnumbers.ntrng
       << setw(8) << gin.randomnumbers.ntgsl
       << "\nEND\n";
  }

  // CENTRE-OF-MASS MOTION

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

  // SPECIAL FORCES

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
  // DIHEDRALRES (promd,md++)
  if(gin.dihedralres.found){
    os << "DIHEDRALRES\n"
       << "#          NTDLR      CDLR    PHILIN\n"
       << setw(16) << gin.dihedralres.ntdlr
       << setw(10) << gin.dihedralres.cdlr
       << setw(10) << gin.dihedralres.philin
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
  // LOCALELEV (promd)
  if(gin.localelev.found){
    os << "LOCALELEV\n"
       << "#    NTLES    NTLEFR    NTLEFU    NLEGRD    NTLESA\n"
       << setw(10) << gin.localelev.ntles
       << setw(10) << gin.localelev.ntlefr
       << setw(10) << gin.localelev.ntlefu
       << setw(10) << gin.localelev.nlegrd
       << setw(10) << gin.localelev.ntlesa
       << "#    CLES       WLES      RLES"
       << setw(10) << gin.localelev.cles
       << setw(10) << gin.localelev.wles
       << setw(10) << gin.localelev.rles
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

  // FREE-ENERGY CALCULATION

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
       << setw(13) << gin.lambdas.ntil
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

  // UMBRELLA (promd)
  if(gin.umbrella.found){
    os << "UMBRELLA\n"
       << "#  NTUS  USCST1  USCST2 USREF1 USREF2\n"
       << setw(7) << gin.umbrella.ntus
       << setw(7) << gin.umbrella.uscst1
       << setw(7) << gin.umbrella.uscst2
       << setw(7) << gin.umbrella.usref1
       << setw(7) << gin.umbrella.usref2
       << "\nEND\n";
  }

  // INPUT-OUTPUT

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
  // EWARN (md++)
  if(gin.ewarn.found){
    os << "EWARN\n"
       << "#  MAXENER\n"
       << setw(10) << gin.ewarn.maxener
       << "\nEND\n";
  }

  // EXTRA

  // Unknown blocks
  for(unsigned int i=0; i< gin.unknown.size(); i++){
    os << gin.unknown[i].name << "\n"
       << gin.unknown[i].content
       << "END\n";
  }


  /*
    Maybe we want to add a string-output to the following blocks?
    LAMBDAS CONSTRAINT MULTIBATH PAIRLIST PRESSURESCALE
  */

    
  return os;
  
}
