// mk_script.h

// We define two global variables. It makes the printing and bookkeeping 
// of the errors and warnings a lot cleaner.
int numWarnings = 0;
int numErrors = 0;
int numTotWarnings = 0;
int numTotErrors = 0;

enum filetype {
  unknownfile, inputfile, topofile, coordfile, refposfile, anatrxfile,
  posresspecfile, xrayfile, disresfile, pttopofile, dihresfile, jvaluefile,
  ledihfile, leumbfile, outputfile, outtrxfile, outtrvfile, outtrffile,
  outtrefile, outtrgfile,
  scriptfile, outbaefile, outbagfile,
  outtrsfile
};
typedef map<string, filetype>::value_type FT;
const FT filetypes[] = {FT("", unknownfile),
  FT("input", inputfile),
  FT("topo", topofile),
  FT("coord", coordfile),
  FT("refpos", refposfile),
  FT("anatrj", anatrxfile),
  FT("posresspec", posresspecfile),
  FT("xray", xrayfile),
  FT("disres", disresfile),
  FT("pttopo", pttopofile),
  FT("dihres", dihresfile),
  FT("jvalue", jvaluefile),
  FT("ledih", ledihfile),
  FT("leumb", leumbfile),
  FT("output", outputfile),
  FT("outtrx", outtrxfile),
  FT("outtrv", outtrvfile),
  FT("outtrf", outtrffile),
  FT("outtre", outtrefile),
  FT("outtrg", outtrgfile),
  FT("outbae", outbaefile),
  FT("outbag", outbagfile),
  FT("outtrs", outtrsfile),
  FT("script", scriptfile)};
const int numFiletypes = sizeof(filetypes)/sizeof(FT);

static map<string, filetype> FILETYPE(filetypes, filetypes + numFiletypes);

enum blocktype {
  unknown, barostatblock, boundcondblock,
  cgrainblock, comtransrotblock, consistencycheckblock,
  constraintblock, covalentformblock, debugblock,
  dihedralresblock, distanceresblock, edsblock, energyminblock,
  ewarnblock, forceblock, geomconstraintsblock,
  gromos96compatblock, initialiseblock, innerloopblock,
  integrateblock, jvalueresblock, lambdasblock,
  localelevblock, multibathblock, multicellblock, multistepblock,
  neighbourlistblock, nonbondedblock, overalltransrotblock,
  pairlistblock, pathintblock, perscaleblock,
  perturbationblock, polariseblock, positionresblock,
  pressurescaleblock, printoutblock, randomnumbersblock,
  readtrajblock, replicablock, rottransblock,
  stepblock, stochdynblock, systemblock,
  thermostatblock, umbrellablock, virialblock,
  writetrajblock, xrayresblock
};

typedef map<string, blocktype>::value_type BT;
int numBlocktypes = 49;
const BT blocktypes[] = {BT("", unknown),
  BT("BAROSTAT", barostatblock),
  BT("BOUNDCOND", boundcondblock),
  BT("CGRAIN", cgrainblock),
  BT("COMTRANSROT", comtransrotblock),
  BT("CONSISTENCYCHECK", consistencycheckblock),
  BT("CONSTRAINT", constraintblock),
  BT("COVALENTFORM", covalentformblock),
  BT("DEBUG", debugblock),
  BT("DIHEDRALRES", dihedralresblock),
  BT("DISTANCERES", distanceresblock),
  BT("EDS", edsblock),
  BT("ENERGYMIN", energyminblock),
  BT("EWARN", ewarnblock),
  BT("FORCE", forceblock),
  BT("GEOMCONSTRAINTS", geomconstraintsblock),
  BT("GROMOS96COMPAT", gromos96compatblock),
  BT("INITIALISE", initialiseblock),
  BT("INNERLOOP", innerloopblock),
  BT("INTEGRATE", integrateblock),
  BT("JVALUERES", jvalueresblock),
  BT("LAMBDAS", lambdasblock),
  BT("LOCALELEV", localelevblock),
  BT("MULTIBATH", multibathblock),
  BT("MULTICELL", multicellblock),
  BT("MULTISTEP", multistepblock),
  BT("NEIGHBOURLIST", neighbourlistblock),
  BT("NONBONDED", nonbondedblock),
  BT("OVERALLTRANSROT", overalltransrotblock),
  BT("PAIRLIST", pairlistblock),
  BT("PATHINT", pathintblock),
  BT("PERSCALE", perscaleblock),
  BT("PERTURBATION", perturbationblock),
  BT("POLARISE", polariseblock),
  BT("POSITIONRES", positionresblock),
  BT("PRESSURESCALE", pressurescaleblock),
  BT("PRINTOUT", printoutblock),
  BT("RANDOMNUMBERS", randomnumbersblock),
  BT("READTRAJ", readtrajblock),
  BT("REPLICA", replicablock),
  BT("ROTTRANS", rottransblock),
  BT("STEP", stepblock),
  BT("STOCHDYN", stochdynblock),
  BT("SYSTEM", systemblock),
  BT("THERMOSTAT", thermostatblock),
  BT("UMBRELLA", umbrellablock),
  BT("VIRIAL", virialblock),
  BT("WRITETRAJ", writetrajblock),
  BT("XRAYRES", xrayresblock)};
static map<string, blocktype> BLOCKTYPE(blocktypes, blocktypes + numBlocktypes);

enum templateelement {
  unknowntemplate, systemtemplate, numbertemplate,
  oldnumbertemplate,
  start_timetemplate, end_timetemplate, queuetemplate
};
typedef map<string, templateelement>::value_type TE;
const TE templateelements[] = {TE("", unknowntemplate),
  TE("system", systemtemplate),
  TE("number", numbertemplate),
  TE("oldnumber", oldnumbertemplate),
  TE("start_time", start_timetemplate),
  TE("end_time", end_timetemplate),
  TE("queue", queuetemplate)};
static map<string, templateelement> TEMPLATE(templateelements,
        templateelements + 7);

//BLOCKDEFINITIONS

class ibarostat {
public:
  int found, ntp, npvar, npcpl[6];
  double comp;

  class pbath {
  public:
    double prsbth;
    vector<double> taubba;
  };
  vector<pbath> pbaths;

  ibarostat() {
    found = 0;
  }
};

class iboundcond {
public:
  int found, ntb, ndfmin;

  iboundcond() {
    found = 0;
  }
};

class icgrain {
public:
  int found, ntcgran;
  double eps;

  icgrain() {
    found = 0;
  }
};

class icomtransrot {
public:
  int found, nscm;

  icomtransrot() {
    found = 0;
  }
};

class iconsistencycheck {
public:
  int found, ntchk, ntckf, ntckv, ntckt;
  int ntcke, ntckr, ntckl;
  double fdckf, fdckv, fdckl;
  vector<int> nckf;

  iconsistencycheck() {
    found = 0;
  }
};

class iconstraint {
public:
  int found, ntc, ntcp, ntcs;
  double ntcp0[3], ntcs0[3];

  iconstraint() {
    found = 0;
    for (int i = 0; i < 0; ++i) {
      ntcp0[i] = ntcs0[i] = -1.0;
    }
  }
};

class icovalentform {
public:
  int found, ntbbh, ntbah, ntbdn;

  icovalentform() {
    found = 0;
  }
};

class idebug {
public:
  int found;

  class routine {
  public:
    string piider;
    int iiideo;
  };
  vector<routine> routines;

  idebug() {
    found = 0;
  }
};

class idihedralres {
public:
  int found, ntdlr;
  double cdlr, philin;

  idihedralres() {
    found = 0;
  }
};

class idistanceres {
public:
  int found, ntdir, ntdira, ntwdir;
  double cdir, dir0, taudir;

  idistanceres() {
    found = 0;
  }
};

class ieds {
public:
  int found, eds, form, numstates;
  vector<double> eir, smooth;
  vector<vector<int> > tree;
  
  ieds() {
    found = 0;
  }
};

class ienergymin {
public:
  int found, ntem, ncyc, nmin;
  double dele, dx0, dxm, flim;

  ienergymin() {
    found = 0;
  }
};

class iewarn {
public:
  int found;
  double maxener;

  iewarn() {
    found = 0;
  }
};

class iforce {
public:
  int found, ntf[10];
  vector<int> nre;

  iforce() {
    found = 0;
  }
};

class igeomconstraints {
public:
  int found, ntcph, ntcpn, ntcs;
  double shktol;

  igeomconstraints() {
    found = 0;
  }
};

class igromos96compat {
public:
  int found, ntnb96, ntr96, ntp96, ntg96;

  igromos96compat() {
    found = 0;
  }
};

class iinitialise {
public:
  int found, ntivel, ntishk, ntinht, ntinhb;
  int ntishi, ntirtc, nticom, ntisti;
  double ig, tempi;

  iinitialise() {
    found = 0;
  }
};

class iinnerloop {
public:
  int found, ntilm, ntils, ntilcd;

  iinnerloop() {
    found = 0;
  }
};

class iintegrate {
public:
  int found, nint;

  iintegrate() {
    found = 0;
  }
};

class ijvalueres {
public:
  int found, ntjvr, ntjvra, le, ngrid, write;
  double cjvr, taujvr, delta;

  ijvalueres() {
    found = 0;
    write = 0;
  }
};

class ilambdas {
public:
  int found, ntil;

  class lambint {
  public:
    int ntli, nilg1, nilg2;
    double ali, bli, cli, dli, eli;
  };
  vector<lambint> lambints;

  ilambdas() {
    found = 0;
  }
};

class ilocalelev {
public:
  int found, ntles, nlepot, ntlesa, ntwle;
  map<int, int> nlepid_ntlerf;

  ilocalelev() {
    found = 0;
  }
};

class imultibath {
public:
  int found, algorithm, num, nbaths, dofset;
  vector<double> temp0, tau;
  vector<int> last, combath, irbath;

  imultibath() {
    found = 0;
    num = -1;
  }
};

class imulticell {
public:
  int found, ntm, ncella, ncellb, ncellc;
  double tolpx, tolpv, tolpf, tolpfw;

  imulticell() {
    found = 0;
  }
};

class imultistep {
public:
  int found, steps, boost;

  imultistep() {
    found = 0;
  }
};

class ineighbourlist {
public:
  int found, plalgo, nupdpl, nupdis, nupdii, type, ncgcen;
  double rcuts, rcuti, gridszx, gridszy, gridszz;

  ineighbourlist() {
    found = 0;
  }
};

class inonbonded {
public:
  int found, nlrele, nshape, na2clc, nkx, nky, nkz, ngx, ngy, ngz;
  int nasord, nfdord, nalias, nqeval, nrdgrd, nwrgrd, nlrlj;
  double appak, rcrf, epsrf, ashape, tola2, epsls, kcut, nspord, faccur, slvdns;

  inonbonded() {
    found = 0;
  }
};

class ioveralltransrot {
public:
  int found, ncmtr, ncmro;
  double cmamx, cmamy, cmamz;

  ioveralltransrot() {
    found = 0;
  }
};

class ipairlist {
public:
  int found, algorithm, nsnb, type;
  double rcutp, rcutl, size;

  ipairlist() {
    found = 0;
  }
};

class ipathint {
public:
  int found, ntpi;

  ipathint() {
    found = 0;
  }
};

class iperscale {
public:
  int found, read, restype;
  double t, kdih, kj, diff, ratio;

  iperscale() {
    found = 0;
  }
};

class iperturbation {
public:
  int found, ntg, nrdgl, nlam, nscale;
  double rlam, dlamt, alphlj, alphc;

  iperturbation() {
    found = 0;
  }
};

class ipolarise {
public:
  int found, cos, efield, damp, write;
  double minfield;

  ipolarise() {
    found = 0;
    write = 0;
  }
};

class ipositionres {
public:
  int found, ntpor, ntporb, ntpors;
  double cpor;

  ipositionres() {
    found = 0;
  }
};

class ipressurescale {
public:
  int found, couple, scale, virial;
  int x_semi, y_semi, z_semi;
  double comp, taup, pres0[3][3];

  ipressurescale() {
    found = 0;
  }
};

class iprintout {
public:
  int found, ntpr, ntpp;

  iprintout() {
    found = 0;
  }
};

class irandomnumbers {
public:
  int found, ntrng, ntgsl;

  irandomnumbers() {
    found = 0;
  }
};

class ireadtraj {
public:
  int found, ntrd, ntrn, ntrb, ntshk;

  ireadtraj() {
    found = 0;
  }
};

class ireplica {
public:
  int found, lrescale, nretrial, nrequil, nrejob, nrewrt;
  vector<double> ret, relam, rets;

  ireplica() {
    found = 0;
    nrewrt = 0;
  }
};

class irottrans {
public:
  int found, rtc, rtclast;

  irottrans() {
    found = 0;
  }
};

class istep {
public:
  int found, nstlim;
  double t, dt;

  istep() {
    found = 0;
  }
};

class istochdyn {
public:
  int found, ntsd, ntfr, nsfr, nbref;
  double rcutf, cfric, tempsd;

  istochdyn() {
    found = 0;
  }
};

class isystem {
public:
  int found, npm, nsm;

  isystem() {
    found = 0;
  }
};

class ithermostat {
public:
  int found, ntt, ntbth, ntset;

  class tbath {
  public:
    int index, ntbtyp, ntbvar;
    double tembth;
    vector<double> taubth;
  };
  vector<tbath> baths;

  class tdofgroup {
  public:
    int ntscpl, ntstyp, ntscns, ntsgt;
    vector<int> ntsgtg;
  };
  vector<tdofgroup> dofgroups;

  ithermostat() {
    found = 0;
  }
};

class iumbrella {
public:
  int found, ntus;
  double uscst1, uscst2, usref1, usref2;

  iumbrella() {
    found = 0;
  }
};

class ivirial {
public:
  int found, ntv, ntvg;

  ivirial() {
    found = 0;
  }
};

class iwritetraj {
public:
  int found, ntwx, ntwse, ntwv, ntwf, ntwe, ntwg, ntwb;

  iwritetraj() {
    found = 0;
    ntwx = 0;
    ntwse = 0;
    ntwv = 0;
    ntwf = 0;
    ntwe = 0;
    ntwg = 0;
    ntwb = 0;
  }
};

class ixrayres {
public:
  int found, ntxr, ntxle, ntwxr, ntwde, ntwxm, rdavg;
  double cxr, cxtau;

  ixrayres() {
    found = 0;
    ntwxr = 0;
  }
};

class iunknown {
public:
  string name;
  string content;

  iunknown(string n) : name(n) {
  }
};

class input {
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
  ieds eds;
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
  imultistep multistep;
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
  ixrayres xrayres;
  vector<iunknown> unknown;
};

class fileInfo {
public:
  gcore::Box box;
  vector<string> blocks;
  vector<int> blockslength;
};

// INSTREAM

// Andreas:
// no checks (except the block size), just read in the data/parameters and check
// if the read in was ok
// more complicated checks are made later, together with the cross checks

istringstream & operator>>(istringstream &is, ibarostat &s) {
  int dum, npbth;
  double taup;
  s.found = 1;
  readValue("BAROSTAT", "NTP", is, s.ntp, "0..3");
  readValue("BAROSTAT", "NPVAR", is, s.npvar, ">=0");
  readValue("BAROSTAT", "NPBTH", is, npbth, ">=0");
  readValue("BAROSTAT", "COMP", is, s.comp, ">=0");
  if (npbth < 0) {
    std::stringstream ss;
    ss << npbth;
    printIO("BAROSTAT", "NPBTH", ss.str(), ">= 0");
  }
  for (int i = 0; i < npbth; ++i) {
    class ibarostat::pbath p;
    is >> dum;
    std::stringstream blockName;
    blockName << "PRSBTH[" << dum << "]";
    readValue("BAROSTAT", blockName.str(), is, p.prsbth, "0,1,2,3");
    if (s.npvar < 0) {
      std::stringstream ss;
      ss << npbth;
      printIO("BAROSTAT", "NPBTH", ss.str(), ">= 0");
    }
    for (int j = 0; j < s.npvar; ++j) {
      std::stringstream blockName;
      blockName << "TAUBBA[" << dum << "," << j + 1 << "]";
      std::stringstream ss;
      ss << taup;
      readValue("BAROSTAT", blockName.str(), is, taup, "> 0.0");
      p.taubba.push_back(taup);
    }
    s.pbaths.push_back(p);
  }
  for (int i = 0; i < 6; i++) {
    std::stringstream blockName;
    blockName << "NPCPL[" << i + 1 << "]";
    readValue("BAROSTAT", blockName.str(), is, s.npcpl[i], "> 0.0");
  }
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of BAROSTAT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iboundcond &s) {
  s.found = 1;
  readValue("BOUNDCOND", "NTB", is, s.ntb, "-1..2");
  readValue("BOUNDCOND", "NDFMIN", is, s.ndfmin, ">= 0");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of BOUNDCOND block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, icgrain &s) {
  s.found = 1;
  readValue("CGRAIN", "NTCGRAN", is, s.ntcgran, "0,1,2");
  readValue("CGRAIN", "EPS", is, s.eps, ">= 0.0");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of CGRAIN block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, icomtransrot &s) {
  s.found = 1;
  readValue("COMTRANSROT", "NSCM", is, s.nscm, "integers");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of COMTRANSROT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iconsistencycheck &s) {
  int nackf, nckf;
  s.found = 1;
  readValue("CONSISTENCYCHECK", "NTCHK", is, s.ntchk, "0, 1");
  readValue("CONSISTENCYCHECK", "NTCKF", is, s.ntckf, "0, 1");
  readValue("CONSISTENCYCHECK", "FDCKF", is, s.fdckf, ">0.0");
  readValue("CONSISTENCYCHECK", "NTCKV", is, s.ntckv, ">0.0");
  readValue("CONSISTENCYCHECK", "FDCKV", is, s.fdckv, ">0.0");
  readValue("CONSISTENCYCHECK", "NTCKT", is, s.ntckt, ">0.0");
  readValue("CONSISTENCYCHECK", "NTCKE", is, s.ntcke, ">0.0");
  readValue("CONSISTENCYCHECK", "NTCKR", is, s.ntckr, ">0.0");
  readValue("CONSISTENCYCHECK", "NTCKL", is, s.ntckl, ">0.0");
  readValue("CONSISTENCYCHECK", "FDCKL", is, s.fdckl, ">0.0");
  readValue("CONSISTENCYCHECK", "NACKF", is, nackf, ">0.0");
  if (nackf < 0) {
    std::stringstream ss;
    ss << nackf;
    printIO("CONSISTENCYCHECK", "NACKF", ss.str(), ">= 0");
  }
  for (int i = 0; i < nackf; i++) {
    readValue("CONSISTENCYCHECK", "NCKF", is, nckf, ">= 1");
    s.nckf.push_back(nckf);
  }
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of CONSISTENCYCHECK block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iconstraint &s) {
  s.found = 1;
  readValue("CONSTRAINT", "NTC", is, s.ntc, "0..4");
  readValue("CONSTRAINT", "NTCP", is, s.ntcp, "1..3");
  readValue("CONSTRAINT", "NTCP0(1)", is, s.ntcp0[0], ">=0");
  if(s.ntcp == 3) {
    readValue("CONSTRAINT", "NTCP0(2)", is, s.ntcp0[1], ">=0");
    readValue("CONSTRAINT", "NTCP0(3)", is, s.ntcp0[2], ">=0");
  }
  readValue("CONSTRAINT", "NTCS", is, s.ntcs, "1..4");
  readValue("CONSTRAINT", "NTCS0(1)", is, s.ntcs0[0], ">=0");
  if (s.ntcs == 3) {
    readValue("CONSTRAINT", "NTCS0(2)", is, s.ntcs0[1], ">=0");
    readValue("CONSTRAINT", "NTCS0(3)", is, s.ntcs0[2], ">=0");
  }
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of CONSTRAINT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, icovalentform &s) {
  s.found = 1;
  readValue("COVALENTFORM", "NTBBH", is, s.ntbbh, "0,1");
  readValue("COVALENTFORM", "NTBAH", is, s.ntbah, "0,1");
  readValue("COVALENTFORM", "NTBDN", is, s.ntbdn, "0,1");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of COVALENTFORM block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, idebug &s) {
  int nrd;
  s.found = 1;
  readValue("DEBUG", "NRDEO", is, nrd, ">=0");
  if (nrd < 0) {
    std::stringstream ss;
    ss << nrd;
    printIO("DEBUG", "NRD", ss.str(), ">= 0");
  }
  for (int i = 0; i < nrd; ++i) {
    class idebug::routine r;
    stringstream ss;
    if (!(is >> r.piider)) {
      stringstream msg;
      msg << "PIIDER(" << i + 1 << ")";
      printIO("DEBUG", msg.str(), r.piider, " a string");
    }
    ss.clear();
    ss.str("");
    stringstream blockName;
    blockName << "IIIDEO(" << i + 1 << ")";
    readValue("DEBUG", blockName.str(), is, r.iiideo, ">=0");
    s.routines.push_back(r);
  }
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of DEBUG block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, idihedralres &s) {
  s.found = 1;
  readValue("DIHEDRALS", "NTDLR", is, s.ntdlr, "0..3");
  readValue("DIHEDRALS", "CDLR", is, s.cdlr, ">=0.0");
  readValue("DIHEDRALS", "PHILIN", is, s.philin, "-1..1");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of DIHEDRALS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, idistanceres &s) {
  s.found = 1;
  readValue("DISTANCERES", "NTDIR", is, s.ntdir, "-2..3");
  readValue("DISTANCERES", "NTDIRA", is, s.ntdira, "0,1");
  readValue("DISTANCERES", "CDIR", is, s.cdir, ">=0.0");
  readValue("DISTANCERES", "DIR0", is, s.dir0, ">=0.0");
  readValue("DISTANCERES", "TAUDIR", is, s.taudir, ">=0.0");
  readValue("DISTANCERES", "NTWDIR", is, s.ntwdir, ">=0");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of DISTANCERES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ieds &s) {
  s.found = 1;
  readValue("EDS", "EDS", is, s.eds, "0,1");
  readValue("EDS", "FORM", is, s.form, "1..3");
  readValue("EDS", "NUMSTATES", is, s.numstates, ">1");
  if (s.numstates <= 1) {
    std::stringstream ss;
    ss << s.numstates;
    printIO("EDS", "NUMSTATES", ss.str(), ">1");
  }
  switch (s.form) {
    case 1:
      s.smooth.resize(1);
      readValue("EDS", "S", is, s.smooth[0], ">0.0");
      break;
    case 2:
      s.smooth.resize(s.numstates*(s.numstates-1)/2);
      for(int i = 0; i < s.numstates*(s.numstates-1)/2; i++) {
        stringstream blockName;
        blockName << "S[" << i + 1 << "]";
        readValue("EDS", blockName.str(), is, s.smooth[i], ">0.0");
      }
      break;
    case 3:
      s.smooth.resize(s.numstates-1);
      s.tree.resize(s.numstates-1);
      for(int N = 0; N < s.numstates-1; N++) {
        stringstream blockName;
        blockName << "i[" << N + 1 << "]";
        vector<int> pair(2);
        readValue("EDS", blockName.str(), is, pair[0], ">0");
        blockName.str("");
        blockName.clear();
        blockName << "j[" << N + 1 << "]";
        readValue("EDS", blockName.str(), is, pair[1], ">0");
        s.tree[N] = pair;
        blockName.str("");
        blockName.clear();
        blockName << "S[" << N + 1 << "]";
        readValue("EDS", blockName.str(), is, s.smooth[N], ">0.0");
      }
      break;
    default:
      stringstream ss;
      ss >> s.numstates;
      printIO("EDS", "NUMSTATES", ss.str(), ">1");
      break;
  }
  s.eir.resize(s.numstates);
  for (int N = 0; N < s.numstates; N++) {
    stringstream blockName;
    blockName << "EIR[" << N + 1 << "]";
    readValue("EDS", blockName.str(), is, s.eir[N], ">0.0");
  }
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of EDS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ienergymin &s) {
  s.found = 1;
  readValue("ENERGYMIN", "NTEM", is, s.ntem, "0..2");
  readValue("ENERGYMIN", "NCYC", is, s.ncyc, ">0");
  readValue("ENERGYMIN", "DELE", is, s.dele, ">0.0");
  readValue("ENERGYMIN", "DX0", is, s.dx0, ">0.0");
  readValue("ENERGYMIN", "DXM", is, s.dxm, ">0.0");
  readValue("ENERGYMIN", "NMIN", is, s.nmin, ">0");
  readValue("ENERGYMIN", "FLIM", is, s.flim, ">0.0");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of ENERGYMIN block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iewarn &s) {
  s.found = 1;
  readValue("EWARN", "MAXENER", is, s.maxener, "a double");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of EWARN block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iforce &s) {
  s.found = 1;
  int negr, nre;
  for (int i = 0; i < 10; i++) {
    readValue("FORCE", "NTF", is, s.ntf[i], "0,1");
  }
  readValue("FORCE", "NEGR", is, negr, ">=0");
  if (negr < 0) {
    std::stringstream ss;
    ss << negr;
    printIO("FORCE", "NEGR", ss.str(), ">=0");
  }
  for (int i = 0; i < negr; i++) {
    stringstream blockName;
    blockName << "NRE(" << i + 1 << ")";
    readValue("FORCE", blockName.str(), is, nre, ">0");
    s.nre.push_back(nre);
  }
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of FORCE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, igeomconstraints &s) {
  s.found = 1;
  readValue("GEOMCONSTRAINTS", "NTCPH", is, s.ntcph, "0,1");
  readValue("GEOMCONSTRAINTS", "NTCPN", is, s.ntcpn, "0,1");
  readValue("GEOMCONSTRAINTS", "NTCS", is, s.ntcs, "0,1");
  readValue("GEOMCONSTRAINTS", "SHKTOL", is, s.shktol, ">0.0");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of GEOMCONSTRAINTS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, igromos96compat &s) {
  s.found = 1;
  readValue("GROMOS96COMPAT", "NTNB96", is, s.ntnb96, "0,1");
  readValue("GROMOS96COMPAT", "NTR96", is, s.ntr96, "0,1");
  readValue("GROMOS96COMPAT", "NTP96", is, s.ntp96, "0,1");
  readValue("GROMOS96COMPAT", "NTG96", is, s.ntg96, "0,1");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of GROMOS96COMPAT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iinitialise &s) {
  s.found = 1;
  readValue("INITIALISE", "NTIVEL", is, s.ntivel, "0,1");
  readValue("INITIALISE", "NTISHK", is, s.ntishk, "0..3");
  readValue("INITIALISE", "NTINHT", is, s.ntinht, "0,1");
  readValue("INITIALISE", "NTINHB", is, s.ntinhb, "0,1");
  readValue("INITIALISE", "NTISHI", is, s.ntishi, "0,1");
  readValue("INITIALISE", "NTIRTC", is, s.ntirtc, "0,1");
  readValue("INITIALISE", "NTICOM", is, s.nticom, "0..3");
  readValue("INITIALISE", "NTISTI", is, s.ntisti, ">0");
  readValue("INITIALISE", "IG", is, s.ig, "0,1");
  readValue("INITIALISE", "TEMPI", is, s.tempi, ">0.0");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of INITIALISE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iinnerloop &s) {
  s.found = 1;
  readValue("INNERLOOP", "NTILM", is, s.ntilm, "0..3");
  readValue("INNERLOOP", "NTILS", is, s.ntils, "0,1");
  readValue("INNERLOOP", "NTILS", is, s.ntilcd, ">=0");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of INNERLOOP block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iintegrate &s) {
  s.found = 1;
  readValue("INTEGRATE", "NINT", is, s.nint, "0,1");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of INTEGRATE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ijvalueres &s) {
  s.found = 1;
  readValue("JVALRES", "NTJVR", is, s.ntjvr, "-3..2");
  readValue("JVALRES", "NTJVRA", is, s.ntjvra, "0..1");
  readValue("JVALRES", "CJVR", is, s.cjvr, ">=0");
  readValue("JVALRES", "TAUJVR", is, s.taujvr, ">=0");
  readValue("JVALRES", "LE", is, s.le, "0..1");
  readValue("JVALRES", "NGRID", is, s.ngrid, ">0");
  readValue("JVALRES", "DELTA", is, s.delta, ">0.0");
  readValue("JVALRES", "NTWJV", is, s.write, ">=0");
  string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of JVALRES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ilambdas &s) {
  s.found = 1;
  readValue("JVALRES", "NTIL", is, s.ntil, "0,1");
  int i = 0;
  string dum;
  while ((is >> dum)) {
    istringstream ss;
    ss.str(dum);
    i++;
    class ilambdas::lambint l;
    stringstream blockName;
    blockName << "NTLI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), ss, l.ntli, "1..11");
    blockName.str("");
    blockName << "NILG1[" << i << "]";
    readValue("LAMBDAS", blockName.str(), ss, l.nilg1, ">0");
    blockName.str("");
    blockName << "NILG2[" << i << "]";
    readValue("LAMBDAS", blockName.str(), ss, l.nilg2, ">0");
    blockName.str("");
    blockName << "ALI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), ss, l.ali, "a double");
    blockName.str("");
    blockName << "BLI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), ss, l.bli, "a double");
    blockName.str("");
    blockName << "CLI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), ss, l.cli, "a double");
    blockName.str("");
    blockName << "DLI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), ss, l.dli, "a double");
    blockName.str("");
    blockName << "ELI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), ss, l.eli, "a double");
  }
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of LAMBDAS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ilocalelev &s) {
  s.found = 1;
  readValue("LOCALELEV", "NTLES", is, s.ntles, "0..2");
  readValue("LOCALELEV", "NLEPOT", is, s.nlepot, ">=0");
  readValue("LOCALELEV", "NTLESA", is, s.ntlesa, "0..2");
  readValue("LOCALELEV", "NTWLE", is, s.ntwle, ">=0");
  if (s.nlepot < 0) {
    std::stringstream ss;
    ss << s.nlepot;
    printIO("LOCALELEV", "NLEPOT", ss.str(), ">=0");
  }
  int nlepid, nlepft;
  string s_nlepid, s_nlepft;
  for (int i = 0; i < s.nlepot; i++) {
    stringstream blockName;
    blockName << "NLEPID(" << i + 1 << ")";
    readValue("LOCALELEV", blockName.str(), is, nlepid, "1..NLEPOT");
    blockName.str("");
    blockName << "NTLEFR(" << i + 1 << ")";
    readValue("LOCALELEV", blockName.str(), is, nlepft, "0,1");
    s.nlepid_ntlerf.insert(pair<int, int> (nlepid, nlepft));
  }
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of LOCALELEV block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, imultibath &s) {
  s.found = 1;
  readValue("MULTIBATH", "ARGORITHM", is, s.algorithm, "0..2");
  if(s.algorithm == 2) {
    readValue("MULTIBATH", "NUM", is, s.num, ">=0");
  }
  readValue("MULTIBATH", "NBATHS", is, s.nbaths, ">=0");
  if (s.nbaths < 0) {
    stringstream ss;
    ss << s.nbaths;
    printIO("MULTIBATH", "NBATHS", ss.str(), ">0");
  }
  for (int i = 0; i < s.nbaths; ++i) {
    double temp0, tau;
    stringstream blockName;
    blockName << "TEMP[" << i + 1 << "]";
    readValue("MULTIBATH", blockName.str(), is, temp0, ">=0.0");
    blockName.str("");
    blockName << "TAU[" << i + 1 << "]";
    readValue("MULTIBATH", blockName.str(), is, tau, ">=0.0");
    s.temp0.push_back(temp0);
    s.tau.push_back(tau);
  }
  readValue("MULTIBATH", "DOFSET", is, s.dofset, ">=0");
  if (s.dofset < 0) {
    stringstream ss;
    ss << s.dofset;
    printIO("MULTIBATH", "DOFSET", ss.str(), ">=0");
  }
  for (int i = 0; i < s.dofset; ++i) {
    int last, combath, irbath;
    stringstream blockName;
    blockName << "LAST[" << i + 1 << "]";
    readValue("MULTIBATH", blockName.str(), is, last, ">=0");
    blockName.str("");
    blockName << "COM-BATH[" << i + 1 << "]";
    readValue("MULTIBATH", blockName.str(), is, combath, ">=1");
    blockName.str("");
    blockName << "IR-BATH[" << i + 1 << "]";
    readValue("MULTIBATH", blockName.str(), is, irbath, ">=1");
    s.last.push_back(last);
    s.combath.push_back(combath);
    s.irbath.push_back(irbath);
  }
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of MULTIBATH block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, imulticell &s) {
  s.found = 1;
  readValue("MULTICELL", "NTM", is, s.ntm, "0,1");
  readValue("MULTICELL", "NCELLA", is, s.ncella, ">=1");
  readValue("MULTICELL", "NCELLB", is, s.ncellb, ">=1");
  readValue("MULTICELL", "NCELLC", is, s.ncellc, ">=1");
  readValue("MULTICELL", "TOLPX", is, s.tolpx, ">=0.0");
  readValue("MULTICELL", "TOLPV", is, s.tolpv, ">=0.0");
  readValue("MULTICELL", "TOLPF", is, s.tolpf, ">=0.0");
  readValue("MULTICELL", "TOLPFW", is, s.tolpfw, ">=0.0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of MULTICELL block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, imultistep &s) {
  s.found = 1;
  readValue("MULTISTEP", "STEPS", is, s.steps, ">0");
  readValue("MULTISTEP", "BOOST", is, s.boost, "0,1");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of MULTISTEP block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ineighbourlist &s) {
  s.found = 1;
  readValue("NEIGHBOURLIST", "PLALGO", is, s.plalgo, "0..2");
  readValue("NEIGHBOURLIST", "NUPDPL", is, s.nupdpl, ">0");
  readValue("NEIGHBOURLIST", "NUPDIS", is, s.nupdis, ">0");
  readValue("NEIGHBOURLIST", "NUPDII", is, s.nupdii, ">0");
  readValue("NEIGHBOURLIST", "RCUTS", is, s.rcuts, ">=0");
  readValue("NEIGHBOURLIST", "RCUTI", is, s.rcuti, ">=RCUTS");
  readValue("NEIGHBOURLIST", "GRIDSZX", is, s.gridszx, ">=0");
  readValue("NEIGHBOURLIST", "GRIDSZY", is, s.gridszy, ">=0");
  readValue("NEIGHBOURLIST", "GRIDSZZ", is, s.gridszz, ">=0");
  readValue("NEIGHBOURLIST", "TYPE", is, s.type, "0..1");
  readValue("NEIGHBOURLIST", "NCGCEN", is, s.ncgcen, ">=-2");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of NEIGHBOURLIST block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, inonbonded &s) {
  s.found = 1;
  readValue("NONBONDED", "NLRELE", is, s.nlrele, "-4..4");
  readValue("NONBONDED", "APPAK", is, s.appak, ">=0.0");
  readValue("NONBONDED", "RCRF", is, s.rcrf, ">=0.0");
  readValue("NONBONDED", "EPSRF", is, s.epsrf, ">=0.0");
  readValue("NONBONDED", "NSHAPE", is, s.nshape, "-1..10");
  readValue("NONBONDED", "ASHAPE", is, s.ashape, ">0.0");
  readValue("NONBONDED", "NA2CLC", is, s.na2clc, "0..4");
  readValue("NONBONDED", "TOLA2", is, s.tola2, ">0.0");
  readValue("NONBONDED", "EPSLS", is, s.epsls, ">0.0");
  readValue("NONBONDED", "NKX", is, s.nkx, ">0");
  readValue("NONBONDED", "NKY", is, s.nky, ">0");
  readValue("NONBONDED", "NKZ", is, s.nkz, ">0");
  readValue("NONBONDED", "KCUT", is, s.kcut, ">0.0");
  readValue("NONBONDED", "NGX", is, s.ngx, ">0");
  readValue("NONBONDED", "NGY", is, s.ngy, ">0");
  readValue("NONBONDED", "NGZ", is, s.ngz, ">0");
  readValue("NONBONDED", "NASORD", is, s.nasord, "1..5");
  readValue("NONBONDED", "NFDORD", is, s.nfdord, "0..5");
  readValue("NONBONDED", "NALIAS", is, s.nalias, ">0");
  readValue("NONBONDED", "NSPORD", is, s.nspord, ">0");
  readValue("NONBONDED", "NQEVAL", is, s.nqeval, ">=0");
  readValue("NONBONDED", "FACCUR", is, s.faccur, ">0.0");
  readValue("NONBONDED", "NRDGRD", is, s.nrdgrd, "0,1");
  readValue("NONBONDED", "NWRGRD", is, s.nwrgrd, "0,1");
  readValue("NONBONDED", "NLRLJ", is, s.nlrlj, "0,1");
  readValue("NONBONDED", "SLVDNS", is, s.slvdns, ">0.0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of NONBONDED block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ioveralltransrot &s) {
  s.found = 1;
  readValue("OVERALLTRANSROT", "NCMTR", is, s.ncmtr, "0,1");
  readValue("OVERALLTRANSROT", "NCMRO", is, s.ncmro, "0,1");
  readValue("OVERALLTRANSROT", "CMAMX", is, s.cmamx, ">=0.0");
  readValue("OVERALLTRANSROT", "CMAMY", is, s.cmamy, ">=0.0");
  readValue("OVERALLTRANSROT", "CMAMZ", is, s.cmamz, ">=0.0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of OVERALLTRANSROT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ipairlist &s) {
  s.found = 1;
  readValue("PAIRLIST", "ALGORITHM", is, s.algorithm, "0,1");
  readValue("PAIRLIST", "NSNB", is, s.nsnb, ">0");
  readValue("PAIRLIST", "RCUTP", is, s.rcutp, ">0.0");
  readValue("PAIRLIST", "RCUTL", is, s.rcutl, ">0.0");
  readValue("PAIRLIST", "SIZE", is, s.size, ">0.0");
  readValue("PAIRLIST", "TYPE", is, s.type, "0,1");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of PAIRLIST block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ipathint &s) {
  s.found = 1;
  readValue("PATHINT", "NTPI", is, s.ntpi, "0,1");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of PATHINT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iperscale &s) {
  s.found = 1;
  readValue("PERSCALE", "RESTYPE", is, s.restype, "0,1");
  readValue("PERSCALE", "KDIH", is, s.kdih, ">=0.0");
  readValue("PERSCALE", "KJ", is, s.kj, ">=0.0");
  readValue("PERSCALE", "T", is, s.t, ">0.0");
  readValue("PERSCALE", "DIFF", is, s.diff, ">=0.0");
  readValue("PERSCALE", "RATIO", is, s.ratio, ">0.0");
  readValue("PERSCALE", "READ", is, s.read, "0,1");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of PERSCALE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iperturbation &s) {
  s.found = 1;
  readValue("PERTURBATION", "NTG", is, s.ntg, "0,1");
  readValue("PERTURBATION", "NRDGL", is, s.nrdgl, "0,1");
  readValue("PERTURBATION", "RLAM", is, s.rlam, "0.0,..,1.0");
  readValue("PERTURBATION", "DLAMT", is, s.dlamt, ">=0.0");
  readValue("PERTURBATION", "ALPHLJ", is, s.alphlj, ">=0.0");
  readValue("PERTURBATION", "ALPHC", is, s.alphc, ">=0.0");
  readValue("PERTURBATION", "NLAM", is, s.nlam, ">0");
  readValue("PERTURBATION", "NSCALE", is, s.nscale, "0..2");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of PERTURBATION block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ipolarise &s) {
  s.found = 1;
  readValue("POLARISE", "COS", is, s.cos, "0,1");
  readValue("POLARISE", "EFIELD", is, s.efield, "0,1");
  readValue("POLARISE", "MINFIELD", is, s.minfield, ">=0.0");
  readValue("POLARISE", "DAMP", is, s.damp, "0,1");
  readValue("POLARISE", "WRITE", is, s.write, ">=0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of POLARISE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ipositionres &s) {
  s.found = 1;
  readValue("POSITIONRES", "NTPOR", is, s.ntpor, "0..3");
  readValue("POSITIONRES", "NTPORB", is, s.ntporb, "0,1");
  readValue("POSITIONRES", "NTPORS", is, s.ntpors, "0,1");
  readValue("POSITIONRES", "CPOR", is, s.cpor, ">=0.0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of POSITIONRES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ipressurescale &s) {
  s.found = 1;
  readValue("PRESSURESCALE", "COUPLE", is, s.couple, "0..2");
  readValue("PRESSURESCALE", "SCALE", is, s.scale, "0..4");
  readValue("PRESSURESCALE", "COMP", is, s.comp, ">0.0");
  readValue("PRESSURESCALE", "TAUP", is, s.taup, ">0.0");
  readValue("PRESSURESCALE", "VIRIAL", is, s.virial, "0..2");
  readValue("PRESSURESCALE", "X_SEMI", is, s.x_semi, "0..2");
  readValue("PRESSURESCALE", "Y_SEMI", is, s.y_semi, "0..2");
  readValue("PRESSURESCALE", "Z_SEMI", is, s.z_semi, "0..2");
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      stringstream blockName;
      blockName << "PRES0[" << i + 1 << "," << j + 1 << "]";
      readValue("PRESSURESCALE", blockName.str(), is, s.pres0[i][j], "a double");
      blockName.str("");
      blockName.clear();
    }
  }
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of PRESSURESCALE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iprintout &s) {
  s.found = 1;
  readValue("PRINTOUT", "NTPR", is, s.ntpr, ">=0");
  readValue("PRINTOUT", "NTPP", is, s.ntpp, "0,1");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of PRINTOUT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, irandomnumbers &s) {
  s.found = 1;
  readValue("RANDOMNUMBERS", "NTRNG", is, s.ntrng, "0,1");
  readValue("RANDOMNUMBERS", "NTGSL", is, s.ntgsl, ">=-1");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of RANDOMNUMBERS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ireadtraj &s) {
  s.found = 1;
  readValue("READTRAJ", "NTRD", is, s.ntrd, "0,1");
  readValue("READTRAJ", "NTRN", is, s.ntrn, "1..18");
  readValue("READTRAJ", "NTRB", is, s.ntrb, "0,1");
  readValue("READTRAJ", "NTSHK", is, s.ntshk, "0,1");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of READTRAJ block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ireplica &s) {
  s.found = 1;
  int nret;
  readValue("REPLICA", "NRET", is, nret, ">=0");
  if (nret < 0) {
    std::stringstream ss;
    ss << nret;
    printIO("REPLICA", "NRET", ss.str(), ">= 0");
  }
  for (int i = 0; i < nret; ++i) {
    stringstream blockName;
    blockName << "RET(" << i + 1 << ")";
    double ret;
    readValue("REPLICA", blockName.str(), is, ret, ">=0.0");
    s.ret.push_back(ret);
  }
  readValue("REPLICA", "LRESCALE", is, s.lrescale, "0,1");
  int nrelam;
  readValue("REPLICA", "NRELAM", is, nrelam, ">=0");
  if (nrelam < 0) {
    std::stringstream ss;
    ss << nrelam;
    printIO("REPLICA", "NRELAM", ss.str(), ">= 0");
  }
  for (int i = 0; i < nrelam; ++i) {
    stringstream blockName;
    blockName << "RELAM(" << i + 1 << ")";
    double relam;
    readValue("REPLICA", blockName.str(), is, relam, ">=0.0");
    s.relam.push_back(relam);
    blockName.str("");
    blockName.clear();
    blockName << "RETS(" << i + 1 << ")";
    double rets;
    readValue("REPLICA", blockName.str(), is, rets, ">=0.0");
    s.rets.push_back(rets);
  }
  readValue("REPLICA", "NRETRIAL", is, s.nretrial, ">=0");
  readValue("REPLICA", "NREQUIL", is, s.nrequil, ">=0");
  readValue("REPLICA", "NREJOB", is, s.nrejob, ">=0");
  readValue("REPLICA", "NREWRT", is, s.nrewrt, ">=0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of REPLICA block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, irottrans &s) {
  s.found = 1;
  readValue("ROTTRANS", "RTC", is, s.rtc, "0,1");
  readValue("ROTTRANS", "RTCLAST", is, s.rtclast, ">0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of ROTTRANS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, istep &s) {
  s.found = 1;
  readValue("STEP", "NSTLIM", is, s.nstlim, ">=0");
  readValue("STEP", "T", is, s.t, ">=0.0");
  readValue("STEP", "DT", is, s.dt, ">0.0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of STEP block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, istochdyn &s) {
  s.found = 1;
  readValue("STOCHDYN", "NTSD", is, s.ntsd, "0,1");
  readValue("STOCHDYN", "NTFR", is, s.ntfr, "0..3");
  readValue("STOCHDYN", "NSFR", is, s.nsfr, ">0");
  readValue("STOCHDYN", "NBREF", is, s.nbref, ">0");
  readValue("STOCHDYN", "RCUTF", is, s.rcutf, ">=0.0");
  readValue("STOCHDYN", "CFRIC", is, s.cfric, ">=0.0");
  readValue("STOCHDYN", "TEMPSD", is, s.tempsd, ">=0.0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of STOCHDYN block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, isystem &s) {
  s.found = 1;
  readValue("SYSTEM", "NPM", is, s.npm, ">=0");
  readValue("SYSTEM", "NSM", is, s.nsm, ">=0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of SYSTEM block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ithermostat &s) {
  s.found = 1;
  readValue("THERMOSTAT", "NTT", is, s.ntt, "0,1");
  readValue("THERMOSTAT", "NTBTH", is, s.ntbth, ">=0");
  readValue("THERMOSTAT", "NTSET", is, s.ntset, ">=0");
  if (s.ntbth < 0) {
    std::stringstream ss;
    ss << s.ntbth;
    printIO("THERMOSTAT", "NTBTH", ss.str(), ">=0");
  }
  for (int i = 0; i < s.ntbth; ++i) {
    class ithermostat::tbath bath;
    stringstream blockName;
    blockName << "I(" << i + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, bath.index, "1..NTBATH");
    blockName.str("");
    blockName.clear();
    blockName << "NTBTYP(" << i + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, bath.ntbtyp, "0..3");
    blockName.str("");
    blockName.clear();
    blockName << "TEMBTH(" << i + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, bath.tembth, ">=0");
    blockName.str("");
    blockName.clear();
    blockName << "NTBVAR(" << i + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, bath.ntbvar, ">=0 and <=MAX_NTBVAR");
    if (bath.ntbvar < 0) {
      std::stringstream ss;
      ss << bath.ntbvar;
      printIO("THERMOSTAT", blockName.str(), ss.str(), ">=0 and <=MAX_NTBVAR");
    }
    for (int j = 0; j < bath.ntbvar; j++) {
      std::stringstream blockName;
      blockName << "TAUBTH(" << i + 1 << "," << j + 1 << ")";
      double tau;
      readValue("THERMOSTAT", blockName.str(), is, tau, ">=0.0");
      bath.taubth.push_back(tau);
    }
    s.baths.push_back(bath);
  }
  if (s.ntset < 0) {
    std::stringstream ss;
    ss << s.ntbth;
    printIO("THERMOSTAT", "NTSET", ss.str(), ">=0");
  }
  for (int j = 0; j < s.ntset; ++j) {
    class ithermostat::tdofgroup dofgroup;
    stringstream blockName;
    blockName << "NTSCPL(" << j + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, dofgroup.ntscpl, "1..NTSET");
    blockName.str("");
    blockName.clear();
    blockName << "NTSTYP(" << j + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, dofgroup.ntstyp, "0..2");
    blockName.str("");
    blockName.clear();
    blockName << "NTSCNS(" << j + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, dofgroup.ntscns, "0..1");
    blockName.str("");
    blockName.clear();
    blockName << "NTSGT(" << j + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, dofgroup.ntsgt, "-2..0 or >=1");
    
    int upper_k = 0;
    if (dofgroup.ntsgt == -2) {
      upper_k = 2;
    } else if (dofgroup.ntsgt == -1) {
      upper_k = 0;
    } else if (dofgroup.ntsgt == 0) {
      stringstream blockName, VAR;
      blockName << "NTSGT(" << j + 1 << ")";
      VAR << dofgroup.ntsgt;
      printIO("THERMOSTAT", blockName.str(), VAR.str(), "not implemented");
    } else { // NTSGT > 0
      upper_k = dofgroup.ntsgt;
    }

    for (int k = 0; k < upper_k; k++) {
      stringstream blockName;
      blockName << "NTSGTG(" << j + 1 << "," << k + 1 << ")";
      int ntsgtg;
      readValue("THERMOSTAT", blockName.str(), is, ntsgtg, "-2..0 or >=1");
      dofgroup.ntsgtg.push_back(ntsgtg);
    }
    s.dofgroups.push_back(dofgroup);
  }
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of THERMOSTAT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iumbrella &s) {
  s.found = 1;
  readValue("UMBRELLA", "NTUS", is, s.ntus, "0,1");
  readValue("UMBRELLA", "USCST1", is, s.uscst1, ">=0");
  readValue("UMBRELLA", "USCST2", is, s.uscst2, ">=0");
  readValue("UMBRELLA", "USREF1", is, s.usref1, ">=0");
  readValue("UMBRELLA", "USREF2", is, s.usref2, ">=0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of UMBRELLA block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ivirial &s) {
  s.found = 1;
  readValue("VIRIAL", "NTV", is, s.ntv, "0,1");
  readValue("VIRIAL", "NTVG", is, s.ntvg, "0..3");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of VIRIAL block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iwritetraj &s) {
  s.found = 1;
  readValue("WRITETRAJ", "NTWX", is, s.ntwx, "an integer");
  readValue("WRITETRAJ", "NTWSE", is, s.ntwse, ">=0");
  readValue("WRITETRAJ", "NTWV", is, s.ntwv, "an integer");
  readValue("WRITETRAJ", "NTWF", is, s.ntwf, "an integer");
  readValue("WRITETRAJ", "NTWE", is, s.ntwe, ">=0");
  readValue("WRITETRAJ", "NTWG", is, s.ntwg, ">=0");
  readValue("WRITETRAJ", "NTWB", is, s.ntwb, ">=0");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of WRITETRAJ block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, ixrayres &s) {
  s.found = 1;
  readValue("XRAYRES", "NTXR", is, s.ntxr, "-2..3");
  readValue("XRAYRES", "NTXLE", is, s.ntxle, "0,1");
  readValue("XRAYRES", "CXR", is, s.cxr, ">=0.0");
  readValue("XRAYRES", "NTWXR", is, s.ntwxr, ">=0");
  readValue("XRAYRES", "NTWDE", is, s.ntwde, "0..3");
  readValue("XRAYRES", "NTWXM", is, s.ntwxm, ">=0");
  readValue("XRAYRES", "CXTAU", is, s.cxtau, ">=0.0");
  readValue("XRAYRES", "RDAVG", is, s.rdavg, "0,1");
  string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      stringstream ss;
      ss << "unexpected end of XRAYRES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

istringstream & operator>>(istringstream &is, iunknown &s) {
  string e = is.str();
  s.content = e.substr(0, e.find("END"));
  return is;
}

Ginstream & operator>>(Ginstream &is, input &gin) {
  vector<string> buffer;
  while (!is.stream().eof()) {
    is.getblock(buffer);

    if (buffer.size() == 1) { // for example if there is twice an "END" in a row
      stringstream msg;
      msg << buffer[0] << " instead of a block name was ignored.";
      printWarning(msg.str());
      buffer.pop_back();
    } else if (buffer.size() > 1) {
      string bufferstring;
      gio::concatenate(buffer.begin() + 1, buffer.end() - 1, bufferstring);
      istringstream bfstream(bufferstring);
      switch (BLOCKTYPE[buffer[0]]) {
        case barostatblock: bfstream >> gin.barostat;
          break;
        case boundcondblock: bfstream >> gin.boundcond;
          break;
        case cgrainblock: bfstream >> gin.cgrain;
          break;
        case comtransrotblock: bfstream >> gin.comtransrot;
          break;
        case consistencycheckblock: bfstream >> gin.consistencycheck;
          break;
        case constraintblock: bfstream >> gin.constraint;
          break;
        case covalentformblock: bfstream >> gin.covalentform;
          break;
        case debugblock: bfstream >> gin.debug;
          break;
        case dihedralresblock: bfstream >> gin.dihedralres;
          break;
        case distanceresblock: bfstream >> gin.distanceres;
          break;
        case energyminblock: bfstream >> gin.energymin;
          break;
        case edsblock: bfstream >> gin.eds;
          break;
        case ewarnblock: bfstream >> gin.ewarn;
          break;
        case forceblock: bfstream >> gin.force;
          break;
        case geomconstraintsblock: bfstream >> gin.geomconstraints;
          break;
        case gromos96compatblock: bfstream >> gin.gromos96compat;
          break;
        case initialiseblock: bfstream >> gin.initialise;
          break;
        case innerloopblock: bfstream >> gin.innerloop;
          break;
        case integrateblock: bfstream >> gin.integrate;
          break;
        case jvalueresblock: bfstream >> gin.jvalueres;
          break;
        case lambdasblock: bfstream >> gin.lambdas;
          break;
        case localelevblock: bfstream >> gin.localelev;
          break;
        case multibathblock: bfstream >> gin.multibath;
          break;
        case multicellblock: bfstream >> gin.multicell;
          break;
        case multistepblock: bfstream >> gin.multistep;
          break;
        case neighbourlistblock: bfstream >> gin.neighbourlist;
          break;
        case nonbondedblock: bfstream >> gin.nonbonded;
          break;
        case overalltransrotblock: bfstream >> gin.overalltransrot;
          break;
        case pairlistblock: bfstream >> gin.pairlist;
          break;
        case pathintblock: bfstream >> gin.pathint;
          break;
        case perscaleblock: bfstream >> gin.perscale;
          break;
        case perturbationblock: bfstream >> gin.perturbation;
          break;
        case polariseblock: bfstream >> gin.polarise;
          break;
        case positionresblock: bfstream >> gin.positionres;
          break;
        case pressurescaleblock: bfstream >> gin.pressurescale;
          break;
        case printoutblock: bfstream >> gin.printout;
          break;
        case randomnumbersblock: bfstream >> gin.randomnumbers;
          break;
        case readtrajblock: bfstream >> gin.readtraj;

          break;
        case replicablock: bfstream >> gin.replica;
          break;
        case rottransblock: bfstream >> gin.rottrans;
          break;
        case stepblock: bfstream >> gin.step;
          break;
        case stochdynblock: bfstream >> gin.stochdyn;
          break;
        case systemblock: bfstream >> gin.system;
          break;
        case thermostatblock: bfstream >> gin.thermostat;
          break;
        case umbrellablock: bfstream >> gin.umbrella;
          break;
        case virialblock: bfstream >> gin.virial;
          break;
        case writetrajblock: bfstream >> gin.writetraj;
          break;
        case xrayresblock: bfstream >> gin.xrayres;
          break;
        case unknown:
          iunknown newblock(buffer[0]);
          bfstream >> newblock;
          gin.unknown.push_back(newblock);
          stringstream msg;
          msg << "Don't know anything about block " << buffer[0]
                  << ". Just storing data.";
          printWarning(msg.str());
      }
    }
  }
  return is;
}

Ginstream & operator>>(Ginstream &is, fileInfo &s) {

  string e;
  string first;
  vector<string> buffer;
  is.getline(first);

  while (!is.stream().eof()) {
    is.getblock(buffer);
    s.blocks.push_back(first);
    s.blockslength.push_back(buffer.size() - 1);
    is.getline(first);
  }
  return is;
}

// TEMPLATE handling of (output) filenames

class filename {
  vector<string> d_parts;
  double d_time, d_dt;
  int d_start;
  string d_system;
  string d_queue;
  string d_template;

public:
  filename();

  filename(string s, double t, double dt, int start = 1, string q = "") {
    d_system = s;
    d_time = t;
    d_dt = dt;
    d_start = start;
    d_queue = q;
  };

  void setInfo(string s, double t, double dt, int start = 1, string q = "") {
    d_system = s;
    d_time = t;
    d_dt = dt;
    d_start = start;
    d_queue = q;
  };

  void setTemplate(string s);

  string temp() {
    return d_template;
  };
  string name(int number);
};

void filename::setTemplate(string s) {
  d_template = s;
  d_parts.clear();
  string::size_type iter;

  string sub;
  iter = s.find('%');
  while (iter != string::npos) {
    sub = s.substr(0, iter);
    s = s.substr(iter + 1, s.size() - iter - 1);
    iter = s.find('%');
    d_parts.push_back(sub);

  }
  d_parts.push_back(s);
}

string filename::name(int number) {
  ostringstream os;
  for (unsigned int i = 0; i < d_parts.size(); i++) {
    if (i % 2) {
      switch (TEMPLATE[d_parts[i]]) {
        case systemtemplate: os << d_system;
          break;
        case numbertemplate: os << d_start + number;
          break;
        case oldnumbertemplate: os << d_start + number - 1;
          break;
        case start_timetemplate: os << d_time + number*d_dt;
          break;
        case end_timetemplate: os << d_time + (number + 1) * d_dt;
          break;
        case queuetemplate: os << d_queue;
          break;
        case unknowntemplate:
          cout << "Do not know how to handle " << d_parts[i]
                  << " in template. Just printing the words." << endl;
          os << d_parts[i];
          break;
      }
    } else os << d_parts[i];
  }
  return os.str();
}

// Jobinfo

class jobinfo {
public:
  map<string, string> param;
  string dir;
  int prev_id;
};


// Writing out of an input file

ostream & operator<<(ostream &os, input &gin) {
  // MOLECULAR SYSTEM

  // SYSTEM (g96, promd, md++)
  if (gin.system.found)
    os << "SYSTEM\n"
          << "#      NPM      NSM\n"
          << setw(10) << gin.system.npm
          << setw(9) << gin.system.nsm
          << "\nEND\n";

  // METHOD EMPLOYED

  // CONSISTENCYCHECK (promd)
  if (gin.consistencycheck.found) {
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
    for (unsigned int i = 0; i < gin.consistencycheck.nckf.size(); ++i) {
      os << setw(10) << gin.consistencycheck.nckf[i];
    }
    os << "\nEND\n";
  }

  // ENERGYMIN (promd, md++): only write if NTEM != 0
  if (gin.energymin.found && gin.energymin.ntem) {
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
  // EDS
  if (gin.eds.found && gin.eds.eds) {
    os << "EDS\n"
            << "#      EDS\n"
            << setw(10) << gin.eds.eds << endl
            << "#     FORM\n"
            << setw(10) << gin.eds.form << endl
            << "#NUMSTATES\n"
            << setw(10) << gin.eds.numstates << endl;
    switch (gin.eds.form) {
      case 1:
        os << "#        S\n"
                << setw(10) << gin.eds.smooth[0] << endl;
        break;
      case 2:
        os << "# S[1..NUMSTATES-1]\n";
        for (int N = 0; N < gin.eds.numstates*(gin.eds.numstates-1)/2; N++) {
          os << setw(10) << gin.eds.smooth[N];
        }
        os << endl;
        break;
      case 3:
        os << "# i[1..NUMSTATES-1]   j[1..NUMSTATES-1]   S[1..NUMSTATES-1]\n";
        for (int N = 0; N < (gin.eds.numstates - 1); N++) {
          os << setw(10) << gin.eds.tree[N][0]
                  << setw(10) << gin.eds.tree[N][1]
                  << setw(10) << gin.eds.smooth[N] << endl;
        }
        break;
    }
    os << "# EIR[1..NUMSTATES]\n";
    for(int N = 0; N < gin.eds.numstates; N++) {
      os << setw(10) << gin.eds.eir[N];
    }
    os << "\nEND\n";
  }
  // STOCHDYN (promd, md++)
  if (gin.stochdyn.found) {
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
  if (gin.readtraj.found) {
    os << "READTRAJ\n"
            << "#     NTRD      NTRN      NTRB     NTSHK\n"
            << setw(10) << gin.readtraj.ntrd
            << setw(10) << gin.readtraj.ntrn
            << setw(10) << gin.readtraj.ntrb
            << setw(10) << gin.readtraj.ntshk
            << "\nEND\n";
  }
  // STEP (promd, md++, g96)
  if (gin.step.found) {
    os << "STEP\n"
            << "#   NSTLIM         T        DT\n"
            << setw(10) << gin.step.nstlim
            << setw(10) << gin.step.t
            << setw(10) << gin.step.dt
            << "\nEND\n";
  }
  // REPLICA (md++)
  if (gin.replica.found) {
    os << "REPLICA\n"
            << "#     NRET\n"
            << setw(10) << gin.replica.ret.size()
            << "\n#  RET(1 ... NRET)\n";
    for (unsigned int i = 0; i < gin.replica.ret.size(); ++i) {
      os << setw(10) << gin.replica.ret[i];
    }
    os << "\n# LRESCALE\n"
            << gin.replica.lrescale
            << "\n#   NRELAM\n"
            << gin.replica.relam.size()
            << "\n#  RELAM(1 ... NRELAM)\n";
    for (unsigned int i = 0; i < gin.replica.relam.size(); ++i) {
      os << setw(10) << gin.replica.relam[i];
    }
    os << "\n#   RETS(1 ... NRELAM)\n";
    for (unsigned int i = 0; i < gin.replica.rets.size(); ++i) {
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
  if (gin.boundcond.found) {
    os << "BOUNDCOND\n"
            << "#      NTB    NDFMIN\n"
            << setw(10) << gin.boundcond.ntb
            << setw(10) << gin.boundcond.ndfmin
            << "\nEND\n";
  }
  // MULTICELL (promd, md++)
  if (gin.multicell.found) {
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
  if (gin.thermostat.found) {
    os << "THERMOSTAT\n"
            << "#       NTT     NTBTH     NTSET\n"
            << setw(11) << gin.thermostat.ntt
            << setw(10) << gin.thermostat.ntbth
            << setw(10) << gin.thermostat.ntset
            << "\n"
            << "# I = 1 ... NTBTH\n"
            << "#         I NTBTYP(I) TEMBTH(I) NTBVAR(I) TAUBTH(I,1...NTVAR)\n";
    for (unsigned int i = 0; i < gin.thermostat.baths.size(); ++i) {
      os << setw(11) << gin.thermostat.baths[i].index
              << setw(10) << gin.thermostat.baths[i].ntbtyp
              << setw(10) << gin.thermostat.baths[i].tembth
              << setw(10) << gin.thermostat.baths[i].ntbvar;
      for (int j = 0; j < gin.thermostat.baths[i].ntbvar; j++) {
        os << setw(10) << gin.thermostat.baths[i].taubth[j];
      }
      os << "\n";
    }
    os << "# NTSCPL(J) NTSTYP(J) NTSCNS(J)  NTSGT(J) NTSGTG(J,1...NTGT(J))\n";
    for (unsigned int j = 0; j < gin.thermostat.dofgroups.size(); ++j) {
      os << setw(11) << gin.thermostat.dofgroups[j].ntscpl
              << setw(10) << gin.thermostat.dofgroups[j].ntstyp
              << setw(10) << gin.thermostat.dofgroups[j].ntscns
              << setw(10) << gin.thermostat.dofgroups[j].ntsgt;
      int upper_k = 0;
      if (gin.thermostat.dofgroups[j].ntsgt == -2) {
        upper_k = 2;
      } else if (gin.thermostat.dofgroups[j].ntsgt == -1) {
        upper_k = 0;
      } else if (gin.thermostat.dofgroups[j].ntsgt == 0) {
        stringstream msg;
        msg << "NTSGT(" << j + 1 << ")";
        printIO("THERMOSTAT", msg.str(), "0", "not implemented");
      } else { // NTSGT > 0
        upper_k = gin.thermostat.dofgroups[j].ntsgt;
      }
      for (int k = 0; k < upper_k; k++) {
        os << setw(10) << gin.thermostat.dofgroups[j].ntsgtg[k];
      }
      os << "\n";
    }
    os << "END\n";
  }
  // MULTIBATH (md++)
  if (gin.multibath.found) {
    os << "MULTIBATH\n"
            << "# ALGORITHM:\n"
            << "#      weak-coupling:      use weak-coupling scheme\n"
            << "#      nose-hoover:        use Nose Hoover scheme\n"
            << "#      nose-hoover-chains: use Nose Hoover chains scheme\n"
            << "# NUM: number of chains in Nose Hoover chains scheme\n"
            << "#      !! only specify NUM when needed !!\n"
            << "# NBATHS: number of temperature baths to couple to\n";
    if (gin.multibath.algorithm == 2) {
      os << "#          ALGORITHM     NUM\n"
              << setw(20) << gin.multibath.algorithm
              << setw(8) << gin.multibath.num;
    } else {
      os << "#          ALGORITHM\n"
              << setw(20) << gin.multibath.algorithm;
    }
    os << "\n#  NBATHS\n"
            << setw(10) << gin.multibath.nbaths
            << "\n";
    os << "# TEMP0(1 ... NBATHS)  TAU(1 ... NBATHS)\n";
    for (int i = 0; i < gin.multibath.nbaths; ++i) {
      os << setw(10) << gin.multibath.temp0[i]
              << setw(10) << gin.multibath.tau[i] << endl;
    }
    os << "\n";
    os << "#   DOFSET: number of distiguishable sets of d.o.f.\n";
    os << setw(10) << gin.multibath.dofset << "\n";
    os << "# LAST(1 ... DOFSET)  COMBATH(1 ... DOFSET)  IRBATH(1 ... DOFSET)\n";
    for (int i = 0; i < gin.multibath.dofset; ++i) {
      os << setw(10) << gin.multibath.last[i]
              << setw(10) << gin.multibath.combath[i]
              << setw(10) << gin.multibath.irbath[i];
    }
    os << "\nEND\n";
  }
  // BAROSTAT (promd)
  if (gin.barostat.found) {
    os << "BAROSTAT\n"
            << "#      NTP      NPVAR     NPBTH      COMP\n"
            << setw(10) << gin.barostat.ntp
            << setw(10) << gin.barostat.npvar
            << setw(10) << gin.barostat.pbaths.size()
            << setw(10) << gin.barostat.comp
            << "\n"
            << "# I = 1 ... NPBTH\n"
            << "# I  PRSBTH(I)  TAUBBA(I,1...NPVAR)\n";
    for (unsigned int i = 0; i < gin.barostat.pbaths.size(); ++i) {
      os << setw(3) << i + 1
              << setw(11) << gin.barostat.pbaths[i].prsbth;
      for (unsigned int j = 0; j < gin.barostat.pbaths[i].taubba.size(); ++j) {
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
  if (gin.virial.found) {
    os << "VIRIAL\n"
            << "#      NTV       NTVG\n"
            << setw(10) << gin.virial.ntv
            << setw(10) << gin.virial.ntvg
            << "\nEND\n";
  }
  // PRESSURESCALE
  if (gin.pressurescale.found) {
    os << "PRESSURESCALE\n"
            << "# COUPLE   SCALE    COMP    TAUP  VIRIAL\n"
            << setw(8) << gin.pressurescale.couple
            << setw(8) << gin.pressurescale.scale << " "
            << setw(8) << gin.pressurescale.comp << " "
            << setw(8) << gin.pressurescale.taup << " "
            << setw(8) << gin.pressurescale.virial
            << "\n# SEMIANISOTROPIC COUPLINGS(X, Y, Z)\n"
            << setw(8) << gin.pressurescale.x_semi << " "
            << setw(8) << gin.pressurescale.y_semi << " "
            << setw(8) << gin.pressurescale.z_semi
            << "\n# PRES0(1...3,1...3)\n"
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
  if (gin.force.found) {
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
    int countnre = 0;
    for (unsigned int i = 0; i < gin.force.nre.size(); i++) {
      os << setw(9) << gin.force.nre[i];
      countnre++;
      if (countnre % 8 == 0) os << endl;
    }
    os << "\nEND\n";
  }
  // MULTISTEP
  if (gin.multistep.found) {
    os << "MULTISTEP\n"
            << "# STEPS calculate non-bonded every STEPSth step.\n"
            << "# BOOST 0,1\n"
            << "#       0: stored forces of STEPSth step are added every step\n"
            << "#       1: stored forces of STEPSth setp are multiplied by STEPS\n"
            << "#       and added every STEPSth step.\n"
            << "#\n"
            << "#" << setw(14) << "STEPS" << setw(15) << "BOOST\n"
            << setw(15) << gin.multistep.steps
            << setw(15) << gin.multistep.boost << endl
            << "END\n";
  }
  // COVALENTFORM (promd, md++)
  if (gin.covalentform.found) {
    os << "COVALENTFORM\n"
            << "#    NTBBH    NTBAH     NTBDN\n"
            << setw(10) << gin.covalentform.ntbbh
            << setw(10) << gin.covalentform.ntbah
            << setw(10) << gin.covalentform.ntbdn
            << "\nEND\n";
  }
  // GEOMCONSTRAINTS (promd)
  if (gin.geomconstraints.found) {
    os << "GEOMCONSTRAINTS\n"
            << "#    NTCPH     NTCPN      NTCS    SHKTOL\n"
            << setw(10) << gin.geomconstraints.ntcph
            << setw(10) << gin.geomconstraints.ntcpn
            << setw(10) << gin.geomconstraints.ntcs
            << setw(10) << gin.geomconstraints.shktol
            << "\nEND\n";
  }
  // CONSTRAINT (md++)
  if (gin.constraint.found) {
    os << "CONSTRAINT\n"
            << "# NTC\n"
            << setw(5) << gin.constraint.ntc
            << "\n";
    if (gin.constraint.ntcp == 3) {
      os << "#      NTCP  NTCP0(1 ... 3)\n"
              << setw(11) << gin.constraint.ntcp
              << setw(10) << gin.constraint.ntcp0[0]
              << setw(10) << gin.constraint.ntcp0[1]
              << setw(10) << gin.constraint.ntcp0[2];
    } else {
      os << "#      NTCP  NTCP0(1)\n"
              << setw(11) << gin.constraint.ntcp
              << setw(10) << gin.constraint.ntcp0[0];
    }
    os << "\n";
    if (gin.constraint.ntcs == 1 || gin.constraint.ntcs == 2 || gin.constraint.ntcs == 5) {
      os << "#      NTCS  NTCS0(1)\n"
              << setw(11) << gin.constraint.ntcs
              << setw(10) << gin.constraint.ntcs0[0];
    } else if (gin.constraint.ntcs == 3) {
      os << "#      NTCS  NTCS0(1 ... 3)\n"
              << setw(11) << gin.constraint.ntcs
              << setw(10) << gin.constraint.ntcs0[0]
              << setw(10) << gin.constraint.ntcs0[1]
              << setw(10) << gin.constraint.ntcs0[2];
    } else {
      os << "#      NTCS\n"
              << setw(11) << gin.constraint.ntcs;
    }
    os << "\nEND\n";
  }
  // GROMOS96COMPAT (promd)
  if (gin.gromos96compat.found) {
    os << "GROMOS96COMPAT\n"
            << "#   NTNB96    NTR96     NTP96     NTG96\n"
            << setw(10) << gin.gromos96compat.ntnb96
            << setw(10) << gin.gromos96compat.ntr96
            << setw(10) << gin.gromos96compat.ntp96
            << setw(10) << gin.gromos96compat.ntg96
            << "\nEND\n";
  }
  // PATHINT (promd)
  if (gin.pathint.found) {
    os << "PATHINT\n"
            << "#  NTPI\n"
            << setw(7) << gin.pathint.ntpi
            << "\nEND\n";
  }
  // POLARISE 
  if (gin.polarise.found) {
    os << "POLARISE\n"
            << "#       COS    EFIELD       MINFIELD      DAMP     WRITE\n"
            << setw(11) << gin.polarise.cos
            << setw(10) << gin.polarise.efield
            << setw(15) << gin.polarise.minfield
            << setw(10) << gin.polarise.damp
            << setw(10) << gin.polarise.write
            << "\nEND\n";
  }
  // INTEGRATE (md++)
  if (gin.integrate.found) {
    os << "INTEGRATE\n"
            << "#   NINT\n"
            << setw(8) << gin.integrate.nint
            << "\nEND\n";
  }
  // CGRAIN (md++)
  if (gin.cgrain.found) {
    os << "CGRAIN\n"
            << "#  NTCGRAN       EPS\n"
            << setw(10) << gin.cgrain.ntcgran
            << setw(10) << gin.cgrain.eps
            << "\nEND\n";
  }
  // ROTTRANS (md++)
  if (gin.rottrans.found) {
    os << "ROTTRANS\n"
            << "#      RTC   RTCLAST\n"
            << setw(10) << gin.rottrans.rtc
            << setw(10) << gin.rottrans.rtclast
            << "\nEND\n";
  }
  // INNERLOOP (md++)
  if (gin.innerloop.found) {
    os << "INNERLOOP\n"
            << "#     NTILM      NTILS\n"
            << setw(10) << gin.innerloop.ntilm
            << setw(10) << gin.innerloop.ntils
            << setw(10) << gin.innerloop.ntilcd
            << "\nEND\n";
  }

  // PAIRLIST GENERATION

  // NEIGHBOURLIST (promd)
  if (gin.neighbourlist.found) {
    os << "NEIGHBOURLIST\n"
            << "#      ALGO  NUPDPL  NUPDIS  NUPDII\n"
            << setw(11) << gin.neighbourlist.plalgo
            << setw(8) << gin.neighbourlist.nupdpl
            << setw(8) << gin.neighbourlist.nupdis
            << setw(8) << gin.neighbourlist.nupdii << endl
            << "#     RCUTS   RCUTI  GRDSZX  GRDSZY  GRDSZZ\n"
            << setw(11) << gin.neighbourlist.rcuts
            << setw(8) << gin.neighbourlist.rcuti
            << setw(8) << gin.neighbourlist.gridszx
            << setw(8) << gin.neighbourlist.gridszy
            << setw(8) << gin.neighbourlist.gridszz << endl
            << "#      TYPE  NCGCEN\n"
            << setw(11) << gin.neighbourlist.type
            << setw(8) << gin.neighbourlist.ncgcen
            << "\nEND\n";
  }
  // PAIRLIST (md++)
  if (gin.pairlist.found) {
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
  if (gin.nonbonded.found) {
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
  if (gin.initialise.found) {
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
  if (gin.randomnumbers.found) {
    os << "RANDOMNUMBERS\n"
            << "#   NTRNG   NTGSL\n"
            << setw(8) << gin.randomnumbers.ntrng
            << setw(8) << gin.randomnumbers.ntgsl
            << "\nEND\n";
  }

  // CENTRE-OF-MASS MOTION

  // OVERALLTRANSROT (promd)
  if (gin.overalltransrot.found) {
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
  if (gin.comtransrot.found) {
    os << "COMTRANSROT\n"
            << "#     NSCM\n"
            << setw(10) << gin.comtransrot.nscm
            << "\nEND\n";
  }

  // SPECIAL FORCES

  // POSITIONRES (promd, md++)
  if (gin.positionres.found) {
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
  // XRAYRES (promd, md++)
  if (gin.xrayres.found) {
    os << "XRAYRES\n"
            << "#    NTXR   NTXLE   CXR   NTWXR   NTWDE   NTWXM   CXTAU  RDAVG\n"
            << setw(10) << gin.xrayres.ntxr
            << setw(10) << gin.xrayres.ntxle
            << setw(10) << gin.xrayres.cxr
            << setw(10) << gin.xrayres.ntwxr
            << setw(10) << gin.xrayres.ntwde
            << setw(10) << gin.xrayres.ntwxm
            << setw(10) << gin.xrayres.cxtau
            << setw(10) << gin.xrayres.rdavg
            << "\nEND\n";
  }
  // DISTANCERES (promd, md++)
  if (gin.distanceres.found) {
    os << "DISTANCERES\n"
            << "# NTDIR\n"
            << "#   0 : no distance restraining\n"
            << "#   -1,1 : use CDIS\n"
            << "#   -2,2: use W0*CDIS\n"
            << "#   NTDIR < 0 : time averaging\n"
            << "#   NTDIR > 0 : no time averaging\n"
            << "# NTDIRA = 1: read in time averaged distances (for continuation run)\n"
            << "# NTDIRA = 0: don't read them in, recalc from scratch\n"
            << "# NTWDIR >= 0 write every NTWDIRth step dist. restr. information to external file\n"
            << "#     NTDIR  NTDIRA    CDIR    DIR0  TAUDIR  NTWDIR\n"
            << setw(11) << gin.distanceres.ntdir
            << setw(8) << gin.distanceres.ntdira
            << setw(8) << gin.distanceres.cdir
            << setw(8) << gin.distanceres.dir0
            << setw(8) << gin.distanceres.taudir
            << setw(8) << gin.distanceres.ntwdir
            << "\nEND\n";
  }
  // DIHEDRALRES (promd,md++)
  if (gin.dihedralres.found) {
    os << "DIHEDRALRES\n"
            << "#          NTDLR      CDLR    PHILIN\n"
            << setw(16) << gin.dihedralres.ntdlr
            << setw(10) << gin.dihedralres.cdlr
            << setw(10) << gin.dihedralres.philin
            << "\nEND\n";
  }
  // JVALUERES (promd, md++)
  if (gin.jvalueres.found) {
    os << "JVALUERES\n"
            << "#        NTJVR  NTJVRA    CJVR  TAUJVR\n"
            << setw(16) << gin.jvalueres.ntjvr
            << setw(10) << gin.jvalueres.ntjvra
            << setw(10) << gin.jvalueres.cjvr
            << setw(10) << gin.jvalueres.taujvr
            << "\n"
            << "#     LE   NGRID    DELTA     NTWJV\n"
            << setw(10) << gin.jvalueres.le
            << setw(10) << gin.jvalueres.ngrid
            << setw(10) << gin.jvalueres.delta
            << setw(10) << gin.jvalueres.write
            << "\nEND\n";
  }
  // LOCALELEV (promd)
  if (gin.localelev.found) {
    os << "LOCALELEV\n"
            << "#     NTLES  NLEPOT  NTLESA    NTWS\n"
            << setw(11) << gin.localelev.ntles
            << setw(8) << gin.localelev.nlepot
            << setw(8) << gin.localelev.ntlesa
            << setw(8) << gin.localelev.ntwle << endl
            << "#    NLEPID       NTLEFR\n";
    for (map<int, int>::iterator it = gin.localelev.nlepid_ntlerf.begin();
            it != gin.localelev.nlepid_ntlerf.end(); ++it) {
      os << setw(10) << it->first << setw(10) << it->second << endl;
    }
    os << "\nEND\n";
  }
  // PERSCALE (md++)
  if (gin.perscale.found) {
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
  if (gin.perturbation.found) {
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
  if (gin.lambdas.found) {
    os << "LAMBDAS\n"
            << "#       NTIL\n"
            << setw(13) << gin.lambdas.ntil
            << "\n# NTLI(1...)  NILG1  NILG2    ALI    BLI    CLI    DLI    ELI\n";
    for (unsigned int i = 0; i < gin.lambdas.lambints.size(); ++i) {
      os << setw(13) << gin.lambdas.lambints[i].ntli
              << setw(7) << gin.lambdas.lambints[i].nilg1
              << setw(7) << gin.lambdas.lambints[i].nilg2
              << setw(7) << gin.lambdas.lambints[i].ali
              << setw(7) << gin.lambdas.lambints[i].bli
              << setw(7) << gin.lambdas.lambints[i].cli
              << setw(7) << gin.lambdas.lambints[i].dli
              << setw(7) << gin.lambdas.lambints[i].eli
              << "\n";
    }
    os << "END\n";
  }

  // UMBRELLA (promd)
  if (gin.umbrella.found) {
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
  if (gin.printout.found) {
    os << "PRINTOUT\n"
            << "#NTPR: print out energies, etc. every NTPR steps\n"
            << "#NTPP: =1 perform dihedral angle transition monitoring\n"
            << "#     NTPR      NTPP\n"
            << setw(10) << gin.printout.ntpr
            << setw(10) << gin.printout.ntpp
            << "\nEND\n";
  }
  // WRITETRAJ (g96)
  if (gin.writetraj.found) {
    os << "WRITETRAJ\n"
            << "#    NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB\n"
            << setw(9) << gin.writetraj.ntwx
            << setw(10) << gin.writetraj.ntwse
            << setw(10) << gin.writetraj.ntwv
            << setw(10) << gin.writetraj.ntwf
            << setw(10) << gin.writetraj.ntwe
            << setw(10) << gin.writetraj.ntwg
            << setw(10) << gin.writetraj.ntwb
            << "\nEND\n";
  }
  // DEBUG (promd)
  if (gin.debug.found) {
    os << "DEBUG\n"
            << "#    NRDEO"
            << setw(10) << gin.debug.routines.size()
            << "\n#  PIIDER(1...NRDEO)  IIIDEO(1...NRDEO)\n";
    for (unsigned int i = 0; i < gin.debug.routines.size(); ++i) {
      os << setw(25) << gin.debug.routines[i].piider
              << setw(8) << gin.debug.routines[i].iiideo
              << "\n";
    }
    os << "END\n";
  }
  // EWARN (md++)
  if (gin.ewarn.found) {
    os << "EWARN\n"
            << "#  MAXENER\n"
            << setw(10) << gin.ewarn.maxener
            << "\nEND\n";
  }

  // EXTRA

  // Unknown blocks
  for (unsigned int i = 0; i < gin.unknown.size(); i++) {
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
