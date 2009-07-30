// mk_script.h

// We define two global variables. It makes the printing and bookkeeping 
// of the errors and warnings a lot cleaner.
int numWarnings = 0;
int numErrors = 0;

enum filetype {
  unknownfile, inputfile, topofile, coordfile, refposfile, anatrxfile,
  posresspecfile, xrayfile, disresfile, pttopofile, dihresfile, jvaluefile,
  ledihfile, outputfile, outtrxfile, outtrvfile, outtrffile,
  outtrefile, outtrgfile,
  scriptfile, outbaefile, outbagfile,
  outtrsfile
};
int numFiletypes = 23;
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
static map<string, filetype> FILETYPE(filetypes, filetypes + numFiletypes);

enum blocktype {
  unknown, barostatblock, boundcondblock,
  cgrainblock, comtransrotblock, consistencycheckblock,
  constraintblock, covalentformblock, debugblock,
  dihedralresblock, distanceresblock, energyminblock,
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
int numBlocktypes = 48;
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
  int found, ntdir, ntdira;
  double cdir, dir0, taudir;

  idistanceres() {
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
  int found, ntilm, ntils;

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
  int found, ntles, nlepot, ntlesa;
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
  int found, nmprpl, nuprpl, nmtwpl, nutwpl, nuirin, nusrin, nmtwin, ncgcen;
  double rcprpl, grprpl, rstwpl, rltwpl, rctwin;

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
  int found, t, read, restype;
  double kdih, kj, diff, ratio;

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
  int found, ntt;

  class tbath {
  public:
    int ntbtyp;
    double tembth;
    vector<double> taubth;
  };
  vector<tbath> baths;

  class dofgroup {
  public:
    int ntscpl, ntstyp, ntscns, ntsgt;
    vector<int> ntsgtg;
  };
  vector<dofgroup> dofgroups;

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
  string e, st;
  int dum, npbth;
  double taup;
  s.found = 1;
  is >> st;
  stringstream ss(st);
  if (!(ss >> s.ntp)) {
    printIO("BAROSTAT", "NTP", st, "0..3");
  }

  is >> s.npvar >> npbth >> s.comp;

  if (s.ntp != 0 && s.ntp != 1 && s.ntp != 2 && s.ntp != 3) {
    std::stringstream ss;
    ss << s.ntp;
    printIO("BAROSTAT", "NTP", ss.str(), "0,1,2,3");
  }
  if (s.npvar < 0) {
    std::stringstream ss;
    ss << s.npvar;
    printIO("BAROSTAT", "NPVAR", ss.str(), ">= 0");
  }
  if (npbth < 0) {
    std::stringstream ss;
    ss << npbth;
    printIO("BAROSTAT", "NPBTH", ss.str(), ">= 0");
  }
  if (s.comp <= 0.0) {
    std::stringstream ss;
    ss << s.comp;
    printIO("BAROSTAT", "COMP", ss.str(), "> 0.0");
  }

  for (int i = 0; i < npbth; ++i) {
    class ibarostat::pbath p;
    is >> dum >> p.prsbth;
    if (p.prsbth < 0.0) {
      std::stringstream si;
      si << "PRSBTH[" << dum << "]";
      std::stringstream ss;
      ss << p.prsbth;
      printIO("BAROSTAT", si.str(), ss.str(), ">= 0.0");
    }
    for (int j = 0; j < s.npvar; ++j) {
      is >> taup;
      p.taubba.push_back(taup);
      if (taup <= 0) {
        std::stringstream si;
        si << "TAUBBA[" << dum << "," << j + 1 << "]";
        std::stringstream ss;
        ss << taup;
        printIO("BAROSTAT", si.str(), ss.str(), "> 0.0");
      }
    }
    s.pbaths.push_back(p);
  }
  for (int i = 0; i < 6; i++) {
    is >> s.npcpl[i];
    if (s.npcpl[i] != 0 && s.npcpl[i] != 1) {
      std::stringstream si;
      si << "NPCPL[" << i + 1 << "]";
      std::stringstream ss;
      ss << s.npcpl[i];
      printIO("BAROSTAT", si.str(), ss.str(), "0,1");
    }
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of BAROSTAT block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iboundcond &s) {
  string e, ntb, ndfmin;
  stringstream ss;
  s.found = 1;
  is >> ntb;
  ss << ntb;
  if (!(ss >> s.ntb)) {
    printIO("BOUNDCOND", "NTB", ntb, "-1..2");
  }
  ss.clear();
  is >> ndfmin;
  ss << ndfmin;
  if (!(ss >> s.ndfmin)) {
    printIO("BOUNDCOND", "NDFMIN", ndfmin, ">= 0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of BOUNDCOND block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, icgrain &s) {
  string e, ntcgran, eps;
  stringstream ss;
  s.found = 1;
  is >> ntcgran;
  ss << ntcgran;
  if (!(ss >> s.ntcgran)) {
    printIO("CGRAIN", "NTCGRAN", ntcgran, "0..2");
  }
  ss.clear();
  is >> eps;
  ss << eps;
  if (!(ss >> s.eps)) {
    printIO("CGRAIN", "EPS", eps, ">= 0.0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of CGRAIN block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, icomtransrot &s) {
  string e, nscm;
  stringstream ss;
  s.found = 1;
  is >> nscm;
  ss << nscm;
  if (!(ss >> s.nscm)) {
    stringstream msg;
    msg << "Could not convert \"" << nscm << "\" to an integer value in COMTRANSROT block";
    printError(msg.str());
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of COMTRANSROT block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iconsistencycheck &s) {
  string e, s_nackf, s_nckf, ntchk, ntckf, fdckf, ntckv, fdckv, ntckt, ntcke, ntckr, ntckl, fdckl;
  int nackf, nckf;
  stringstream ss;
  s.found = 1;
  is >> ntchk;
  ss << ntchk;
  if (!(ss >> s.ntchk)) {
    printIO("CONSISTENCYCHECK", "NTCHK", ntchk, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> ntckf;
  ss << ntckf;
  if (!(ss >> s.ntckf)) {
    printIO("CONSISTENCYCHECK", "NTCKF", ntckf, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> fdckf;
  ss << fdckf;
  if (!(ss >> s.fdckf)) {
    printIO("CONSISTENCYCHECK", "FDCKF", fdckf, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> ntckv;
  ss << ntckv;
  if (!(ss >> s.ntckv)) {
    printIO("CONSISTENCYCHECK", "NTCKV", ntckv, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> fdckv;
  ss << fdckv;
  if (!(ss >> s.fdckv)) {
    printIO("CONSISTENCYCHECK", "FDCKV", fdckv, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> ntckt;
  ss << ntckt;
  if (!(ss >> s.ntckt)) {
    printIO("CONSISTENCYCHECK", "NTCKT", ntckt, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> ntcke;
  ss << ntcke;
  if (!(ss >> s.ntcke)) {
    printIO("CONSISTENCYCHECK", "NTCKE", ntcke, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> ntckr;
  ss << ntckr;
  if (!(ss >> s.ntckr)) {
    printIO("CONSISTENCYCHECK", "NTCKR", ntckr, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> ntckl;
  ss << ntckl;
  if (!(ss >> s.ntckl)) {
    printIO("CONSISTENCYCHECK", "NTCKL", ntckl, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> fdckl;
  ss << fdckl;
  if (!(ss >> s.fdckl)) {
    printIO("CONSISTENCYCHECK", "FDCKL", fdckl, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> s_nackf;
  ss << s_nackf;
  if (!(ss >> nackf)) {
    printIO("CONSISTENCYCHECK", "NACKF", s_nackf, ">= 0");
  }
  ss.clear();
  ss.str("");
  if (nackf < 0) {
    std::stringstream ss;
    ss << nackf;
    printIO("CONSISTENCYCHECK", "NACKF", ss.str(), ">= 0");
  }
  for (int i = 0; i < nackf; i++) {
    is >> s_nckf;
    ss << s_nckf;
    if (!(ss >> nckf)) {
      printIO("CONSISTENCYCHECK", "NCKF", s_nckf, ">= 1");
    }
    ss.clear();
    ss.str("");
    s.nckf.push_back(nckf);
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of CONSISTENCYCHECK block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iconstraint &s) {
  string e, ntc, ntcp, ntcp0, ntcs, ntcs0;
  stringstream ss;
  s.found = 1;
  is >> ntc;
  ss << ntc;
  if (!(ss >> s.ntc)) {
    printIO("CONSTRAINT", "NTC", ntc, "0..4");
  }
  ss.clear();
  ss.str("");
  is >> ntcp;
  ss << ntcp;
  if (!(ss >> s.ntcp)) {
    printIO("CONSTRAINT", "NTCP", ntcp, "1..3");
  }
  ss.clear();
  ss.str("");
  is >> ntcp0;
  ss << ntcp0;
  if (!(ss >> s.ntcp0[0])) {
    printIO("CONSTRAINT", "NTCP0(1)", ntcp0, ">= 0");
  }
  ss.clear();
  ss.str("");
  if (s.ntcp == 3) {
    is >> ntcp;
    ss << ntcp;
    if (!(ss >> s.ntcp0[1])) {
      printIO("CONSTRAINT", "NTCP0(2)", ntcp, ">= 0");
    }
    ss.clear();
    ss.str("");
    is >> ntcp;
    ss << ntcp;
    if (!(ss >> s.ntcp0[2])) {
      printIO("CONSTRAINT", "NTCP0(3)", ntcp, ">= 0");
    }
    ss.clear();
    ss.str("");
  }
  is >> ntcs;
  ss << ntcs;
  if (!(ss >> s.ntcs)) {
    printIO("CONSTRAINT", "NTCS", ntcs, "1..4");
  }
  ss.clear();
  ss.str("");
  if (s.ntcs != 4) {
    is >> ntcs0;
    ss << ntcs0;
    if (!(ss >> s.ntcs0[0])) {
      printIO("CONSTRAINT", "NTCS0(1)", ntcs0, ">= 0");
    }
    ss.clear();
    ss.str("");
  }
  if (s.ntcs == 3) {
    is >> ntcs0;
    ss << ntcs0;
    if (!(ss >> s.ntcs0[1])) {
      printIO("CONSTRAINT", "NTCS0(2)", ntcs0, ">= 0");
    }
    ss.clear();
    ss.str("");
    is >> ntcs0;
    ss << ntcs0;
    if (!(ss >> s.ntcs0[2])) {
      printIO("CONSTRAINT", "NTCS0(3)", ntcs0, ">= 0");
    }
    ss.clear();
    ss.str("");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of CONSTRAINT block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, icovalentform &s) {
  string e, ntbbh, ntbah, ntbdn;
  stringstream ss;
  s.found = 1;
  is >> ntbbh;
  ss << ntbbh;
  if (!(ss >> s.ntbbh)) {
    printIO("COVALENTFORM", "NTBBH", ntbbh, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> ntbah;
  ss << ntbah;
  if (!(ss >> s.ntbah)) {
    printIO("COVALENTFORM", "NTBAH", ntbah, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> ntbdn;
  ss << ntbdn;
  if (!(ss >> s.ntbdn)) {
    printIO("COVALENTFORM", "NTBDN", ntbdn, "0, 1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of COVALENTFORM block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, idebug &s) {
  string e, s_nrd, piider, iiideo;
  stringstream ss;
  int nrd;
  s.found = 1;
  is >> s_nrd;
  ss << s_nrd;
  if (!(ss >> nrd)) {
    printIO("DEBUG", "NRDEO", s_nrd, ">= 0");
  }
  ss.clear();
  ss.str("");
  if (nrd < 0) {
    std::stringstream ss;
    ss << nrd;
    printIO("DEBUG", "NRD", ss.str(), ">= 0");
  }
  for (int i = 0; i < nrd; ++i) {
    class idebug::routine r;
    is >> piider;
    ss << piider;
    if (!(ss >> r.piider)) {
      stringstream msg;
      msg << "PIIDER(" << i + 1 << ")";
      printIO("DEBUG", msg.str(), piider, " a string");
    }
    ss.clear();
    ss.str("");
    is >> iiideo;
    ss << iiideo;
    if (!(ss >> r.iiideo)) {
      stringstream msg;
      msg << "IIIDEO(" << i + 1 << ")";
      printIO("DEBUG", msg.str(), iiideo, ">= 0");
    }
    ss.clear();
    ss.str("");
    s.routines.push_back(r);
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of DEBUG block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, idihedralres &s) {
  string e, ntdlr, cdlr, philin;
  stringstream ss;
  s.found = 1;
  is >> ntdlr;
  ss << ntdlr;
  if (!(ss >> s.ntdlr)) {
    printIO("DIHEDRALS", "NTDLR", ntdlr, "0..3");
  }
  ss.clear();
  ss.str("");
  is >> cdlr;
  ss << cdlr;
  if (!(ss >> s.cdlr)) {
    printIO("DIHEDRALS", "CDLR", cdlr, ">= 0.0");
  }
  ss.clear();
  ss.str("");
  is >> philin;
  ss << philin;
  if (!(ss >> s.philin)) {
    printIO("DIHEDRALS", "PHILIN", philin, "-1..1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of DIHEDRALRES block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, idistanceres &s) {
  string e, ntdir, ntdira, cdir, dir0, taudir;
  stringstream ss;
  s.found = 1;
  is >> ntdir;
  ss << ntdir;
  if (!(ss >> s.ntdir)) {
    printIO("DISTANCERES", "NTDIR", ntdir, "-2..3");
  }
  ss.clear();
  ss.str("");
  is >> ntdira;
  ss << ntdira;
  if (!(ss >> s.ntdira)) {
    printIO("DISTANCERES", "NTDIRA", ntdira, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> cdir;
  ss << cdir;
  if (!(ss >> s.cdir)) {
    printIO("DISTANCERES", "CDIR", cdir, ">= 0.0");
  }
  ss.clear();
  ss.str("");
  is >> dir0;
  ss << dir0;
  if (!(ss >> s.dir0)) {
    printIO("DISTANCERES", "DIR0", dir0, ">= 0.0");
  }
  ss.clear();
  ss.str("");
  is >> taudir;
  ss << taudir;
  if (!(ss >> s.taudir)) {
    printIO("DISTANCERES", "TAUDIR", taudir, ">= 0.0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of DISTANCERES block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ienergymin &s) {
  string e, ntem, ncyc, dele, dx0, dxm, nmin, flim;
  stringstream ss;
  s.found = 1;
  is >> ntem;
  ss << ntem;
  if (!(ss >> s.ntem)) {
    printIO("ENERGYMIN", "NTEM", ntem, "0..2");
  }
  ss.clear();
  ss.str("");
  is >> ncyc;
  ss << ncyc;
  if (!(ss >> s.ncyc)) {
    printIO("ENERGYMIN", "NCYC", ncyc, "> 0");
  }
  ss.clear();
  ss.str("");
  is >> dele;
  ss << dele;
  if (!(ss >> s.dele)) {
    printIO("ENERGYMIN", "DELE", dele, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> dx0;
  ss << dx0;
  if (!(ss >> s.dx0)) {
    printIO("ENERGYMIN", "DX0", dx0, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> dxm;
  ss << dxm;
  if (!(ss >> s.dxm)) {
    printIO("ENERGYMIN", "DXM", dxm, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> nmin;
  ss << nmin;
  if (!(ss >> s.nmin)) {
    printIO("ENERGYMIN", "NMIN", nmin, "> 0");
  }
  ss.clear();
  ss.str("");
  is >> flim;
  ss << flim;
  if (!(ss >> s.flim)) {
    printIO("ENERGYMIN", "FLIM", flim, "> 0.0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of ENERGYMIN block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iewarn &s) {
  string e, maxener;
  stringstream ss;
  s.found = 1;
  is >> maxener;
  ss << maxener;
  if (!(ss >> s.maxener)) {
    printIO("EWARN", "MAXENER", maxener, "a double");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of EWARN block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, iforce &s) {
  string e, ntf, s_negr, s_nre;
  stringstream ss;
  s.found = 1;
  int negr, nre;
  for (int i = 0; i < 10; i++) {
    is >> ntf;
    ss << ntf;
    if (!(ss >> s.ntf[i])) {
      printIO("FORCE", "NTF", ntf, "0, 1");
    }
    ss.clear();
    ss.str("");
  }
  is >> s_negr;
  ss << s_negr;
  if (!(ss >> negr)) {
    printIO("FORCE", "NEGR", s_negr, ">= 0");
  }
  ss.clear();
  ss.str("");
  if (negr <= 0) {
    std::stringstream ss;
    ss << negr;
    printIO("FORCE", "NEGR", ss.str(), "> 0");
  }
  for (int i = 0; i < negr; i++) {
    is >> s_nre;
    ss << s_nre;
    if (!(ss >> nre)) {
      stringstream msg;
      msg << "NRE(" << i + 1 << ")";
      printIO("FORCE", msg.str(), s_nre, "> 0");
    }
    ss.clear();
    ss.str("");
    s.nre.push_back(nre);
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of FORCE block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, igeomconstraints &s) {
  string e, ntcph, ntcpn, ntcs, shktol;
  stringstream ss;
  s.found = 1;
  is >> ntcph;
  ss << ntcph;
  if (!(ss >> s.ntcph)) {
    printIO("GEOMCONSTRAINTS", "NTCPH", ntcph, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> ntcpn;
  ss << ntcpn;
  if (!(ss >> s.ntcpn)) {
    printIO("GEOMCONSTRAINTS", "NTCPN", ntcpn, "0, 1");
  }
  ss.clear();
  ss.str("");
  if (!(ss >> s.ntcs)) {
    printIO("GEOMCONSTRAINTS", "NTCS", ntcs, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> shktol;
  ss << shktol;
  if (!(ss >> s.shktol)) {
    printIO("GEOMCONSTRAINTS", "SHKTOL", shktol, "> 0.0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of GEOMCONSTRAINTS block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, igromos96compat &s) {
  string e, ntnb96, ntr96, ntp96, ntg96;
  stringstream ss;
  s.found = 1;
  is >> ntnb96;
  ss << ntnb96;
  if (!(ss >> s.ntnb96)) {
    printIO("GROMOS96COMPAT", "NTNB96", ntnb96, "0, 1");
    cout << "ss = " << ss.str() << endl;
  }
  ss.clear();
  ss.str("");
  ss.str("");
  is >> ntr96;
  ss << ntr96;
  if (!(ss >> s.ntr96)) {
    printIO("GROMOS96COMPAT", "NTR96", ntr96, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> ntp96;
  ss << ntp96;
  if (!(ss >> s.ntp96)) {
    printIO("GROMOS96COMPAT", "NTP96", ntp96, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> ntg96;
  ss << ntg96;
  if (!(ss >> s.ntg96)) {
    printIO("GROMOS96COMPAT", "NTG96", ntg96, "0, 1");
  }
  is >> e;

  if (e != "") {
    stringstream ss;
    ss << "unexpected end of GROMOS96COMPAT block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, iinitialise &s) {
  string e, ntivel, ntishk, ntinht, ntinhb, ntishi, ntirtc, nticom, ntisti, ig, tempi;
  stringstream ss;
  s.found = 1;
  is >> ntivel;
  ss << ntivel;
  if (!(ss >> s.ntivel)) {
    printIO("INITIALISE", "NTIVEL", ntivel, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> ntishk;
  ss << ntishk;
  if (!(ss >> s.ntishk)) {
    printIO("INITIALISE", "NTISHK", ntishk, "0..3");
  }
  ss.clear();
  ss.str("");
  is >> ntinht;
  ss << ntinht;
  if (!(ss >> s.ntinht)) {
    printIO("INITIALISE", "NTINHT", ntinht, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> ntinhb;
  ss << ntinhb;
  if (!(ss >> s.ntinhb)) {
    printIO("INITIALISE", "NTINHB", ntinhb, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> ntishi;
  ss << ntishi;
  if (!(ss >> s.ntishi)) {
    printIO("INITIALISE", "NTISHI", ntishi, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> ntirtc;
  ss << ntirtc;
  if (!(ss >> s.ntirtc)) {
    printIO("INITIALISE", "NTIRTC", ntirtc, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> nticom;
  ss << nticom;
  if (!(ss >> s.nticom)) {
    printIO("INITIALISE", "NTICOM", nticom, "0..3");
  }
  ss.clear();
  ss.str("");
  is >> ntisti;
  ss << ntisti;
  if (!(ss >> s.ntisti)) {
    printIO("INITIALISE", "NTISTI", ntisti, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> ig;
  ss << ig;
  if (!(ss >> s.ig)) {
    printIO("INITIALISE", "IG", ig, "> 0");
  }
  ss.clear();
  ss.str("");
  is >> tempi;
  ss << tempi;
  if (!(ss >> s.tempi)) {
    printIO("INITIALISE", "TEMPI", tempi, "> 0.0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of INITIALISE block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iinnerloop &s) {
  string e, ntilm, ntils;
  stringstream ss;
  s.found = 1;
  is >> ntilm;
  ss << ntilm;
  if (!(ss >> s.ntilm)) {
    printIO("INNERLOOP", "NTILM", ntilm, "0..3");
  }
  ss.clear();
  ss.str("");
  is >> ntils;
  ss << ntils;
  if (!(ss >> s.ntils)) {
    printIO("INNERLOOP", "NTILS", ntils, "0, 1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of INNERLOOP block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iintegrate &s) {
  string e, nint;
  stringstream ss;
  s.found = 1;
  is >> nint;
  ss << nint;
  if (!(ss >> s.nint)) {
    printIO("INTEGRATE", "NINT", nint, "0, 1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of INTEGRATE block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ijvalueres &s) {
  string e, ntjvr, ntjvra, cjvr, taujvr, le, ngrid, delta, write;
  stringstream ss;
  s.found = 1;
  is >> ntjvr;
  ss << ntjvr;
  if (!(ss >> s.ntjvr)) {
    printIO("JVALRES", "NTJVR", ntjvr, "-3..2");
  }
  ss.clear();
  ss.str("");
  is >> ntjvra;
  ss << ntjvra;
  if (!(ss >> s.ntjvra)) {
    printIO("JVALRES", "NTJVRA", ntjvra, "0..1");
  }
  ss.clear();
  ss.str("");
  is >> cjvr;
  ss << cjvr;
  if (!(ss >> s.cjvr)) {
    printIO("JVALRES", "CJVR", cjvr, ">= 0");
  }
  ss.clear();
  ss.str("");
  is >> taujvr;
  ss << taujvr;
  if (!(ss >> s.taujvr)) {
    printIO("JVALRES", "TAUJVR", taujvr, ">= 0");
  }
  ss.clear();
  ss.str("");
  is >> le;
  ss << le;
  if (!(ss >> s.le)) {
    printIO("JVALRES", "LE", le, "0..1");
  }
  ss.clear();
  ss.str("");
  is >> ngrid;
  ss << ngrid;
  if (!(ss >> s.ngrid)) {
    printIO("JVALRES", "NGRID", ngrid, "> 0");
  }
  ss.clear();
  ss.str("");
  is >> delta;
  ss << delta;
  if (!(ss >> s.delta)) {
    printIO("JVALRES", "DELTA", delta, "> 0.0");
  }
  ss.clear();
  ss.str("");
  is >> write;
  ss << write;
  if (!(ss >> s.write)) {
    printIO("JVALRES", "NTWJV", write, ">= 0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of JVALUERES block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ilambdas &s) {
  string e, dum;
  s.found = 1;
  is >> s.ntil;
  if (s.ntil != 0 && s.ntil != 1) {
    std::stringstream ss;
    ss << s.ntil;
    printIO("LAMBDAS", "NTIL", ss.str(), "off(0),on(1)");
  }
  int i = 0;
  while ((is >> dum)) {
    i++;
    std::transform(dum.begin(), dum.end(), dum.begin(), ::tolower);
    class ilambdas::lambint l;
    if (dum == "bond") l.ntli = 1;
    else if (dum == "angle") l.ntli = 2;
    else if (dum == "dihedral") l.ntli = 3;
    else if (dum == "improper") l.ntli = 4;
    else if (dum == "vdw") l.ntli = 5;
    else if (dum == "vdw_soft") l.ntli = 6;
    else if (dum == "crf") l.ntli = 7;
    else if (dum == "crf_soft") l.ntli = 8;
    else if (dum == "distanceres") l.ntli = 9;
    else if (dum == "dihedralres") l.ntli = 10;
    else if (dum == "mass") l.ntli = 11;
    else {
      std::stringstream ss(dum);
      if ((!(ss >> l.ntli)) || l.ntli < 1 || l.ntli > 11) {
        std::stringstream si;
        si << "NTLI[" << i << "]";
        printIO("LAMBDAS", si.str(), dum,
                "bond(1), angle(2), dihedral(3), improper(4), vdw(5), "
                "vdw_soft(6), crf(7), crf_soft(8), distanceres(9), "
                "dihedralres(10), mass(11)");
      }
    }
    is >> l.nilg1 >> l.nilg2 >> l.ali >> l.bli >> l.cli >> l.dli >> l.eli;
    if (l.nilg1 <= 0) {
      std::stringstream ss, si;
      ss << l.nilg1;
      si << "NILG1[" << i << "]";
      printIO("LAMBDAS", si.str(), ss.str(), "> 0");
    }
    if (l.nilg2 <= 0) {
      std::stringstream ss, si;
      ss << l.nilg2;
      si << "NILG2[" << i << "]";
      printIO("LAMBDAS", si.str(), ss.str(), "> 0");
    }

    s.lambints.push_back(l);

  }
  e = dum;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of LAMBDAS block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ilocalelev &s) {
  string e, ntles, nlepot, ntlesa;
  stringstream ss;
  s.found = 1;
  is >> ntles;
  ss << ntles;
  if (!(ss >> s.ntles)) {
    printIO("LOCALELEV", "NTLES", ntles, "0..2");
  }
  ss.clear();
  ss.str("");
  is >> nlepot;
  ss << nlepot;
  if (!(ss >> s.nlepot)) {
    printIO("LOCALELEV", "NLEPOT", nlepot, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> ntlesa;
  ss << ntlesa;
  if (!(ss >> s.ntlesa)) {
    printIO("LOCALELEV", "NTLESA", ntlesa, "0..2");
  }
  ss.clear();
  ss.str("");
  if (s.nlepot < 0) {
    std::stringstream ss;
    ss << s.nlepot;
    printIO("LOCALELEV", "NLEPOT", ss.str(), ">= 0");
  }
  int id, erfl;
  string s_id, s_erfl;
  for (int i = 0; i < s.nlepot; i++) {
    is >> s_id;
    ss << s_id;
    if (!(ss >> id)) {
      stringstream msg;
      msg << "NLEPID(" << i + 1 << ")";
      printIO("LOCALELEV", msg.str(), s_id, "1..NLEPOT");
    }
    ss.clear();
    ss.str("");
    is >> s_erfl;
    ss << s_erfl;
    if (!(ss >> erfl)) {
      stringstream msg;
      msg << "NTLEFR(" << i + 1 << ")";
      printIO("LOCALELEV", msg.str(), s_erfl, "0, 1");
    }
    ss.clear();
    ss.str("");
    s.nlepid_ntlerf.insert(pair<int, int> (id, erfl));
  }
  if (int(s.nlepid_ntlerf.size()) != s.nlepot) {
    std::stringstream ss, si;
    ss << s.nlepid_ntlerf.size() << " potentials";
    si << s.nlepot;
    printIO("LOCALELEV", "NLEPID and NTLEFR", ss.str(), si.str());
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of LOCALELEV block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, imultibath &s) {
  string e, algorithm, num, nbaths, s_temp0, s_tau, dofset, s_last, s_combath, s_irbath;
  stringstream ss;
  double temp0, tau;
  int last, combath, irbath;
  s.found = 1;
  is >> algorithm;
  ss << algorithm;
  if (!(ss >> s.algorithm)) {
    printIO("MULTIBATH", "ARGORITHM", algorithm, "0..2");
  }
  ss.clear();
  ss.str("");
  if (s.algorithm == 2) {
    is >> num;
    ss << num;
    if (!(ss >> s.num)) {
      printIO("MULTIBATH", "NUM", algorithm, ">= 0");
    }
    ss.clear();
    ss.str("");
  }
  is >> nbaths;
  ss << nbaths;
  if (!(ss >> s.nbaths)) {
    printIO("MULTIBATH", "NBATHS", nbaths, ">= 0");
  }
  ss.clear();
  ss.str("");
  if (s.nbaths < 0) {
    std::stringstream ss;
    ss << s.nbaths;
    printIO("MULTIBATH", "NBATHS", ss.str(), "> 0");
  }
  for (int i = 0; i < s.nbaths; ++i) {
    is >> s_temp0;
    ss << s_temp0;
    if (!(ss >> temp0)) {
      stringstream msg;
      msg << "TEMP(" << i + 1 << ")";
      printIO("MULTIBATH", msg.str(), s_temp0, ">= 0.0");
    }
    ss.clear();
    ss.str("");
    is >> s_tau;
    ss << s_tau;
    if (!(ss >> tau)) {
      stringstream msg;
      msg << "TAU(" << i + 1 << ")";
      printIO("MULTIBATH", msg.str(), s_tau, ">= 0.0");
    }
    ss.clear();
    ss.str("");
    s.temp0.push_back(temp0);
    s.tau.push_back(tau);
  }
  is >> dofset;
  ss << dofset;
  if (!(ss >> s.dofset)) {
    printIO("MULTIBATH", "DOFSET", dofset, ">= 0");
  }
  ss.clear();
  ss.str("");
  if (s.dofset < 0) {
    stringstream ss;
    ss << s.dofset;
    printIO("MULTIBATH", "DOFSET", ss.str(), ">= 0");
  }
  for (int i = 0; i < s.dofset; ++i) {
    is >> s_last;
    ss << s_last;
    if (!(ss >> last)) {
      stringstream msg;
      msg << "LAST(" << i + 1 << ")";
      printIO("MULTIBATH", msg.str(), s_last, ">= 0");
    }
    ss.clear();
    ss.str("");
    is >> s_combath;
    ss << s_combath;
    if (!(ss >> combath)) {
      stringstream msg;
      msg << "COM-BATH(" << i + 1 << ")";
      printIO("MULTIBATH", msg.str(), s_combath, ">= 1");
    }
    ss.clear();
    ss.str("");
    ss.clear();
    ss.str("");
    is >> s_irbath;
    ss << s_irbath;
    if (!(ss >> irbath)) {
      stringstream msg;
      msg << "IR-BATH(" << i + 1 << ")";
      printIO("MULTIBATH", msg.str(), s_irbath, ">= 1");
    }
    ss.clear();
    ss.str("");
    s.last.push_back(last);
    s.combath.push_back(combath);
    s.irbath.push_back(irbath);
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of MULTIBATH block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, imulticell &s) {
  string e, ntm, ncella, ncellb, ncellc, tolpx, tolpv, tolpf, tolpfw;
  stringstream ss;
  s.found = 1;
  is >> ntm;
  ss << ntm;
  if (!(ss >> s.ntm)) {
    printIO("MULTICELL", "NTM", ntm, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> ncella;
  ss << ncella;
  if (!(ss >> s.ncella)) {
    printIO("MULTICELL", "NCELLA", ncella, ">=1");
  }
  ss.clear();
  ss.str("");
  is >> ncellb;
  ss << ncellb;
  if (!(ss >> s.ncellb)) {
    printIO("MULTICELL", "NCELLB", ncellb, ">=1");
  }
  ss.clear();
  ss.str("");
  is >> ncellc;
  ss << ncellc;
  if (!(ss >> s.ncellc)) {
    printIO("MULTICELL", "NCELLC", ncellc, ">=1");
  }
  ss.clear();
  ss.str("");
  is >> tolpx;
  ss << tolpx;
  if (!(ss >> s.tolpx)) {
    printIO("MULTICELL", "TOLPX", tolpx, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> tolpv;
  ss << tolpv;
  if (!(ss >> s.tolpv)) {
    printIO("MULTICELL", "TOLPV", tolpv, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> tolpf;
  ss << tolpf;
  if (!(ss >> s.tolpf)) {
    printIO("MULTICELL", "TOLPF", tolpf, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> tolpfw;
  ss << tolpfw;
  if (!(ss >> s.tolpfw)) {
    printIO("MULTICELL", "TOLPFW", tolpfw, ">0.0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of MULTICELL block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, imultistep &s) {
  string e, steps, boost;
  stringstream ss;
  s.found = 1;
  is >> steps;
  ss << steps;
  if (!(ss >> s.steps)) {
    printIO("MULTISTEP", "STEPS", steps, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> boost;
  ss << boost;
  if (!(ss >> s.boost)) {
    printIO("MULTISTEP", "BOOST", boost, "0,1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of MULTISTEP block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

// Andreas:
// promd block, still needs to be changed but wait for after the GROMOS meeting
// on July 2 2009
istringstream & operator>>(istringstream &is, ineighbourlist &s) {
  string e, nmprpl, nuprpl, rcprpl, grprpl, nmtwpl, nutwpl, rstwpl;
  string rltwpl, nuirin, nusrin, nmtwin, rctwin, ncgcen;
  stringstream ss;
  s.found = 1;
  is >> nmprpl;
  ss << nmprpl;
  if (!(ss >> s.nmprpl)) {
    printIO("NEIGHBOURLIST", "NMPRPL", nmprpl, "");
  }
  ss.clear();
  ss.str("");
  is >> nuprpl;
  ss << nuprpl;
  if (!(ss >> s.nuprpl)) {
    printIO("NEIGHBOURLIST", "NUPRPL", nuprpl, "");
  }
  ss.clear();
  ss.str("");
  is >> rcprpl;
  ss << rcprpl;
  if (!(ss >> s.rcprpl)) {
    printIO("NEIGHBOURLIST", "RCPRPL", rcprpl, "");
  }
  ss.clear();
  ss.str("");
  is >> grprpl;
  ss << grprpl;
  if (!(ss >> s.grprpl)) {
    printIO("NEIGHBOURLIST", "GRPRPL", grprpl, "");
  }
  ss.clear();
  ss.str("");
  is >> nmtwpl;
  ss << nmtwpl;
  if (!(ss >> s.nmtwpl)) {
    printIO("NEIGHBOURLIST", "NMTWPL", nmtwpl, "");
  }
  ss.clear();
  ss.str("");
  is >> nutwpl;
  ss << nutwpl;
  if (!(ss >> s.nutwpl)) {
    printIO("NEIGHBOURLIST", "NUTWPL", nutwpl, "");
  }
  ss.clear();
  ss.str("");
  is >> rstwpl;
  ss << rstwpl;
  if (!(ss >> s.rstwpl)) {
    printIO("NEIGHBOURLIST", "RSTWPL", rstwpl, "");
  }
  ss.clear();
  ss.str("");
  is >> rltwpl;
  ss << rltwpl;
  if (!(ss >> s.rltwpl)) {
    printIO("NEIGHBOURLIST", "RLTWPL", rltwpl, "");
  }
  ss.clear();
  ss.str("");
  is >> nuirin;
  ss << nuirin;
  if (!(ss >> s.nuirin)) {
    printIO("NEIGHBOURLIST", "NUIRIN", nuirin, "");
  }
  ss.clear();
  ss.str("");
  is >> nusrin;
  ss << nusrin;
  if (!(ss >> s.nusrin)) {
    printIO("NEIGHBOURLIST", "NUSRIN", nusrin, "");
  }
  ss.clear();
  ss.str("");
  is >> nmtwin;
  ss << nmtwin;
  if (!(ss >> s.nmtwin)) {
    printIO("NEIGHBOURLIST", "NMTWIN", nmtwin, "");
  }
  ss.clear();
  ss.str("");
  is >> rctwin;
  ss << rctwin;
  if (!(ss >> s.rctwin)) {
    printIO("NEIGHBOURLIST", "RCTWIN", rctwin, "");
  }
  ss.clear();
  ss.str("");
  is >> ncgcen;
  ss << ncgcen;
  if (!(ss >> s.ncgcen)) {
    printIO("NEIGHBOURLIST", "NCGCEN", ncgcen, ">= -2");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of NEIGBOURLIST block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, inonbonded &s) {
  string e, nlrele, appak, rcrf, epsrf, nshape, ashape, na2clc;
  string tola2, epsls, nkx, nky, nkz, kcut, ngx, ngy, ngz, nasord;
  string nfdord, nalias, nspord, nqeval, faccur, nrdgrd, nwrgrd, nlrlj, slvdns;
  stringstream ss;
  s.found = 1;
  is >> nlrele;
  ss << nlrele;
  if (!(ss >> s.nlrele)) {
    printIO("NONBONDED", "NLRELE", nlrele, "-4..4");
  }
  ss.clear();
  ss.str("");
  is >> appak;
  ss << appak;
  if (!(ss >> s.appak)) {
    printIO("NONBONDED", "APPAK", appak, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> rcrf;
  ss << rcrf;
  if (!(ss >> s.rcrf)) {
    printIO("NONBONDED", "RCRF", rcrf, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> epsrf;
  ss << epsrf;
  if (!(ss >> s.epsrf)) {
    printIO("NONBONDED", "EPSRF", epsrf, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> nshape;
  ss << nshape;
  if (!(ss >> s.nshape)) {
    printIO("NONBONDED", "NSHAPE", nshape, "-1..10");
  }
  ss.clear();
  ss.str("");
  is >> ashape;
  ss << ashape;
  if (!(ss >> s.ashape)) {
    printIO("NONBONDED", "ASHAPE", ashape, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> na2clc;
  ss << na2clc;
  if (!(ss >> s.na2clc)) {
    printIO("NONBONDED", "NA2CLC", na2clc, "0..4");
  }
  ss.clear();
  ss.str("");
  is >> tola2;
  ss << tola2;
  if (!(ss >> s.tola2)) {
    printIO("NONBONDED", "TOLA2", tola2, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> epsls;
  ss << epsls;
  if (!(ss >> s.epsls)) {
    printIO("NONBONDED", "EPSLS", epsls, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> nkx;
  ss << nkx;
  if (!(ss >> s.nkx)) {
    printIO("NONBONDED", "NKX", nkx, ">0");
  }
  ss.clear();
  ss.str("");
  is >> nky;
  ss << nky;
  if (!(ss >> s.nky)) {
    printIO("NONBONDED", "NKY", nky, ">0");
  }
  ss.clear();
  ss.str("");
  is >> nkz;
  ss << nkz;
  if (!(ss >> s.nkz)) {
    printIO("NONBONDED", "NKZ", nkz, ">0");
  }
  ss.clear();
  ss.str("");
  is >> kcut;
  ss << kcut;
  if (!(ss >> s.kcut)) {
    printIO("NONBONDED", "KCUT", kcut, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> ngx;
  ss << ngx;
  if (!(ss >> s.ngx)) {
    printIO("NONBONDED", "NGX", ngx, ">0");
  }
  ss.clear();
  ss.str("");
  is >> ngy;
  ss << ngy;
  if (!(ss >> s.ngy)) {
    printIO("NONBONDED", "NGY", ngy, ">0");
  }
  ss.clear();
  ss.str("");
  is >> ngz;
  ss << ngz;
  if (!(ss >> s.ngz)) {
    printIO("NONBONDED", "NGZ", ngz, ">0");
  }
  ss.clear();
  ss.str("");
  is >> nasord;
  ss << nasord;
  if (!(ss >> s.nasord)) {
    printIO("NONBONDED", "NASORD", nasord, "1..5");
  }
  ss.clear();
  ss.str("");
  is >> nfdord;
  ss << nfdord;
  if (!(ss >> s.nfdord)) {
    printIO("NONBONDED", "NFDORD", nfdord, "0..5");
  }
  ss.clear();
  ss.str("");
  is >> nalias;
  ss << nalias;
  if (!(ss >> s.nalias)) {
    printIO("NONBONDED", "NALIAS", nalias, ">0");
  }
  ss.clear();
  ss.str("");
  is >> nspord;
  ss << nspord;
  if (!(ss >> s.nspord)) {
    printIO("NONBONDED", "NSPORD", nspord, ">0");
  }
  ss.clear();
  ss.str("");
  is >> nqeval;
  ss << nqeval;
  if (!(ss >> s.nqeval)) {
    printIO("NONBONDED", "NQEVAL", nqeval, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> faccur;
  ss << faccur;
  if (!(ss >> s.faccur)) {
    printIO("NONBONDED", "FACCUR", faccur, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> nrdgrd;
  ss << nrdgrd;
  if (!(ss >> s.nrdgrd)) {
    printIO("NONBONDED", "NRDGRD", nrdgrd, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> nwrgrd;
  ss << nwrgrd;
  if (!(ss >> s.nwrgrd)) {
    printIO("NONBONDED", "NWRGRD", nwrgrd, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> nlrlj;
  ss << nlrlj;
  if (!(ss >> s.nlrlj)) {
    printIO("NONBONDED", "NLRLJ", nlrlj, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> slvdns;
  ss << slvdns;
  if (!(ss >> s.slvdns)) {
    printIO("NONBONDED", "SLVDNS", slvdns, ">0.0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of NONBONDED block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ioveralltransrot &s) {
  string e, ncmtr, ncmro, cmamx, cmamy, cmamz;
  stringstream ss;
  s.found = 1;
  is >> ncmtr;
  ss << ncmtr;
  if (!(ss >> s.ncmtr)) {
    printIO("OVERALLTRANSROT", "NCMTR", ncmtr, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> ncmro;
  ss << ncmro;
  if (!(ss >> s.ncmro)) {
    printIO("OVERALLTRANSROT", "NCMRO", ncmro, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> cmamx;
  ss << cmamx;
  if (!(ss >> s.cmamx)) {
    printIO("OVERALLTRANSROT", "CMAMX", cmamx, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> cmamy;
  ss << cmamy;
  if (!(ss >> s.cmamy)) {
    printIO("OVERALLTRANSROT", "CMAMY", cmamy, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> cmamz;
  ss << cmamz;
  if (!(ss >> s.cmamz)) {
    printIO("OVERALLTRANSROT", "CMAMZ", cmamz, ">=0.0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of OVERALLTRANSROT block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ipairlist &s) {
  string e, algorithm, nsnb, rcutp, rcutl, size, type;
  stringstream ss;
  s.found = 1;
  is >> algorithm;
  ss << algorithm;
  if (!(ss >> s.algorithm)) {
    printIO("PAIRLIST", "ALGORITHM", algorithm, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> nsnb;
  ss << nsnb;
  if (!(ss >> s.nsnb)) {
    printIO("PAIRLIST", "NSNB", nsnb, ">0");
  }
  ss.clear();
  ss.str("");
  is >> rcutp;
  ss << rcutp;
  if (!(ss >> s.rcutp)) {
    printIO("PAIRLIST", "RCUTP", rcutp, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> rcutl;
  ss << rcutl;
  if (!(ss >> s.rcutl)) {
    printIO("PAIRLIST", "RCUTL", rcutl, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> size;
  ss << size;
  if (!(ss >> s.size)) {
    printIO("PAIRLIST", "SIZE", size, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> type;
  ss << type;
  if (!(ss >> s.type)) {
    printIO("PAIRLIST", "TYPE", type, "0,1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of PAIRLIST block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ipathint &s) {
  string e, ntpi;
  stringstream ss;
  s.found = 1;
  is >> ntpi;
  ss << ntpi;
  if (!(ss >> s.ntpi)) {
    printIO("PATHINT", "NTPI", ntpi, "0,1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of PATHINT block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iperscale &s) {
  string e, restype, kdih, kj, t, diff, ratio, read;
  stringstream ss;
  s.found = 1;
  is >> restype;
  ss << restype;
  if (!(ss >> s.restype)) {
    printIO("PERSCALE", "RESTYPE", restype, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> kdih;
  ss << kdih;
  if (!(ss >> s.kdih)) {
    printIO("PERSCALE", "KDIH", kdih, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> kj;
  ss << kj;
  if (!(ss >> s.kj)) {
    printIO("PERSCALE", "KJ", kj, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> t;
  ss << t;
  if (!(ss >> s.t)) {
    printIO("PERSCALE", "T", t, ">0");
  }
  ss.clear();
  ss.str("");
  is >> diff;
  ss << diff;
  if (!(ss >> s.diff)) {
    printIO("PERSCALE", "DIFF", diff, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> ratio;
  ss << ratio;
  if (!(ss >> s.ratio)) {
    printIO("PERSCALE", "RATIO", ratio, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> read;
  ss << read;
  if (!(ss >> s.read)) {
    printIO("PERSCALE", "READ", read, "0,1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of PERSCALE block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iperturbation &s) {
  string e, ntg, nrdgl, rlam, dlamt, alphlj, alphc, nlam, nscale;
  stringstream ss;
  s.found = 1;
  is >> ntg;
  ss << ntg;
  if (!(ss >> s.ntg)) {
    printIO("PERTURBATION", "NTG", ntg, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> nrdgl;
  ss << nrdgl;
  if (!(ss >> s.nrdgl)) {
    printIO("PERTURBATION", "NRDGL", nrdgl, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> rlam;
  ss << rlam;
  if (!(ss >> s.rlam)) {
    printIO("PERTURBATION", "RLAM", rlam, "0.0..1.0");
  }
  ss.clear();
  ss.str("");
  is >> dlamt;
  ss << dlamt;
  if (!(ss >> s.dlamt)) {
    printIO("PERTURBATION", "DLAMT", dlamt, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> alphlj;
  ss << alphlj;
  if (!(ss >> s.alphlj)) {
    printIO("PERTURBATION", "ALPHLJ", alphlj, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> alphc;
  ss << alphc;
  if (!(ss >> s.alphc)) {
    printIO("PERTURBATION", "ALPHC", alphc, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> nlam;
  ss << nlam;
  if (!(ss >> s.nlam)) {
    printIO("PERTURBATION", "NLAM", nlam, ">0");
  }
  ss.clear();
  ss.str("");
  is >> nscale;
  ss << nscale;
  if (!(ss >> s.nscale)) {
    printIO("PERTURBATION", "NSCALE", nscale, "0..2");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of PERTURBATION block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ipolarise &s) {
  string e, cos, efield, minfield, damp, write;
  stringstream ss;
  s.found = 1;
  is >> cos;
  ss << cos;
  if (!(ss >> s.cos)) {
    printIO("POLARISE", "COS", cos, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> efield;
  ss << efield;
  if (!(ss >> s.efield)) {
    printIO("POLARISE", "EFIELD", efield, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> minfield;
  ss << minfield;
  if (!(ss >> s.minfield)) {
    printIO("POLARISE", "MINFIELD", minfield, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> damp;
  ss << damp;
  if (!(ss >> s.damp)) {
    printIO("POLARISE", "DAMP", damp, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> write;
  ss << write;
  if (!(ss >> s.write)) {
    printIO("POLARISE", "WRITE", write, ">=0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of POLARISE block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ipositionres &s) {
  string e, ntpor, ntporb, ntpors, cpor;
  stringstream ss;
  s.found = 1;
  is >> ntpor;
  ss << ntpor;
  if (!(ss >> s.ntpor)) {
    printIO("POSITIONRES", "NTPOR", ntpor, "0..3");
  }
  ss.clear();
  ss.str("");
  is >> ntporb;
  ss << ntporb;
  if (!(ss >> s.ntporb)) {
    printIO("POSITIONRES", "NTPORB", ntporb, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> ntpors;
  ss << ntpors;
  if (!(ss >> s.ntpors)) {
    printIO("POSITIONRES", "NTPORS", ntpors, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> cpor;
  ss << cpor;
  if (!(ss >> s.cpor)) {
    printIO("POSITIONRES", "CPOR", cpor, ">=0.0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of POSITIONRES block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ipressurescale &s) {
  string e, couple, scale, comp, taup, virial;
  string pres0;
  stringstream ss;
  s.found = 1;
  is >> couple;
  ss << couple;
  if (!(ss >> s.couple)) {
    printIO("PRESSURESCALE", "COUPLE", couple, "0..2");
  }
  ss.clear();
  ss.str("");
  is >> scale;
  ss << scale;
  if (!(ss >> s.scale)) {
    printIO("PRESSURESCALE", "SCALE", scale, "0..3");
  }
  ss.clear();
  ss.str("");
  is >> comp;
  ss << comp;
  if (!(ss >> s.comp)) {
    printIO("PRESSURESCALE", "COMP", comp, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> taup;
  ss << taup;
  if (!(ss >> s.taup)) {
    printIO("PRESSURESCALE", "TAUP", taup, ">0.0");
  }
  ss.clear();
  ss.str("");
  is >> virial;
  ss << virial;
  if (!(ss >> s.virial)) {
    printIO("PRESSURESCALE", "VIRIAL", virial, "0..2");
  }
  ss.clear();
  ss.str("");
  is >> pres0;
  ss << pres0;
  if(!(ss >> s.pres0[0][0])){
    printIO("PRESSURESCALE", "PRES0(1,1)", pres0, ": a numerical value");
  }
  ss.clear();
  ss.str("");
  is >> pres0;
  ss << pres0;
  if(!(ss >> s.pres0[0][1])){
    printIO("PRESSURESCALE", "PRES0(1,2)", pres0, ": a numerical value");
  }
  ss.clear();
  ss.str("");
  is >> pres0;
  ss << pres0;
  if(!(ss >> s.pres0[0][2])){
    printIO("PRESSURESCALE", "PRES0(1,3)", pres0, ": a numerical value");
  }
  ss.clear();
  ss.str("");
  is >> pres0;
  ss << pres0;
  if(!(ss >> s.pres0[1][0])){
    printIO("PRESSURESCALE", "PRES0(2,1)", pres0, ": a numerical value");
  }
  ss.clear();
  ss.str("");
  is >> pres0;
  ss << pres0;
  if(!(ss >> s.pres0[1][1])){
    printIO("PRESSURESCALE", "PRES0(2,2)", pres0, ": a numerical value");
  }
  ss.clear();
  ss.str("");
  is >> pres0;
  ss << pres0;
  if(!(ss >> s.pres0[1][2])){
    printIO("PRESSURESCALE", "PRES0(2,3)", pres0, ": a numerical value");
  }
  ss.clear();
  ss.str("");
  is >> pres0;
  ss << pres0;
  if(!(ss >> s.pres0[2][0])){
    printIO("PRESSURESCALE", "PRES0(3,1)", pres0, ": a numerical value");
  }
  ss.clear();
  ss.str("");
  is >> pres0;
  ss << pres0;
  if(!(ss >> s.pres0[2][1])){
    printIO("PRESSURESCALE", "PRES0(3,2)", pres0, ": a numerical value");
  }
  ss.clear();
  ss.str("");
  is >> pres0;
  ss << pres0;
  if(!(ss >> s.pres0[2][2])){
    printIO("PRESSURESCALE", "PRES0(3,3)", pres0, ": a numerical value");
  }
  ss.clear();
  ss.str("");
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of PRESSURESCALE block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iprintout &s) {
  string e, ntpr, ntpp;
  stringstream ss;
  s.found = 1;
  is >> ntpr;
  ss << ntpr;
  if (!(ss >> s.ntpr)) {
    printIO("PRINTOUT", "NTPR", ntpr, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> ntpp;
  ss << ntpp;
  if (!(ss >> s.ntpp)) {
    printIO("PRINTOUT", "NTPP", ntpp, "0,1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of PRINTOUT block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, irandomnumbers &s) {
  string e, ntrng, ntgsl;
  stringstream ss;
  s.found = 1;
  is >> ntrng;
  ss << ntrng;
  if (!(ss >> s.ntrng)) {
    printIO("RANDOMNUMBERS", "NTRNG", ntrng, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> ntgsl;
  ss << ntgsl;
  if (!(ss >> s.ntgsl)) {
    printIO("RANDOMNUMBERS", "NTGSL", ntgsl, ">= -1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of RANDOMNUMBERS block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ireadtraj &s) {
  string e, ntrd, ntrn, ntrb, ntshk;
  stringstream ss;
  s.found = 1;
  is >> ntrd;
  ss << ntrd;
  if (!(ss >> s.ntrd)) {
    printIO("READTRAJ", "NTRD", ntrd, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> ntrn;
  ss << ntrn;
  if (!(ss >> s.ntrn)) {
    printIO("READTRAJ", "NTRN", ntrn, "1..18");
  }
  ss.clear();
  ss.str("");
  is >> ntrb;
  ss << ntrb;
  if (!(ss >> s.ntrb)) {
    printIO("READTRAJ", "NTRB", ntrb, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> ntshk;
  ss << ntshk;
  if (!(ss >> s.ntshk)) {
    printIO("READTRAJ", "NTSHK", ntshk, "0,1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of READTRAJ block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ireplica &s) {
  string e, s_nret, s_ret, s_rets, lrescale, s_nrelam, s_relam, nretrial, nrequil, nrejob, nrewrt;
  stringstream ss;
  int nret, nrelam;
  double ret, rets, relam;
  s.found = 1;
  is >> s_nret;
  ss << s_nret;
  if (!(ss >> nret)) {
    printIO("REPLICA", "NRET", s_nret, ">= 0");
  }
  ss.clear();
  ss.str("");
  if (nret < 0) {
    std::stringstream ss;
    ss << nret;
    printIO("REPLICA", "NRET", ss.str(), ">= 0");
  }
  for (int i = 0; i < nret; ++i) {
    is >> ret;
    ss << ret;
    if (!(ss >> ret)) {
      stringstream msg;
      msg << "RET(" << i + 1 << ")";
      printIO("REPLICA", msg.str(), s_ret, ">= 0.0");
    }
    ss.clear();
    ss.str("");
    s.ret.push_back(ret);
  }
  is >> lrescale;
  ss << lrescale;
  if (!(ss >> s.lrescale)) {
    printIO("REPLICA", "LRESCALE", lrescale, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> s_nrelam;
  ss << s_nrelam;
  if (!(ss >> nrelam)) {
    printIO("REPLICA", "NRELAM", s_nrelam, ">= 0.0");
  }
  ss.clear();
  ss.str("");
  if (nrelam < 0) {
    std::stringstream ss;
    ss << nrelam;
    printIO("REPLICA", "NRELAM", ss.str(), ">= 0");
  }
  for (int i = 0; i < nrelam; ++i) {
    is >> s_relam;
    ss << s_relam;
    if (!(ss >> relam)) {
      stringstream msg;
      msg << "RELAM(" << i + 1 << ")";
      printIO("REPLICA", msg.str(), s_relam, ">= 0.0");
    }
    ss.clear();
    ss.str("");
    s.relam.push_back(relam);


    is >> s_rets;
    ss << s_rets;
    if (!(ss >> rets)) {
      stringstream msg;
      msg << "RETS(" << i + 1 << ")";
      printIO("REPLICA", msg.str(), s_rets, ">= 0.0");
    }
    s.rets.push_back(rets);
  }
  ss.clear();
  ss.str("");
  is >> nretrial;
  ss << nretrial;
  if (!(ss >> s.nretrial)) {
    printIO("REPLICA", "NRETRIAL", nretrial, ">= 0");
  }
  ss.clear();
  ss.str("");
  is >> nrequil;
  ss << nrequil;
  if (!(ss >> s.nrequil)) {
    printIO("REPLICA", "NREQUIL", nrequil, ">= 0");
  }
  ss.clear();
  ss.str("");
  is >> nrejob;
  ss << nrejob;
  if (!(ss >> s.nrejob)) {
    printIO("REPLICA", "NREJOB", nrejob, "> 0");
  }
  ss.clear();
  ss.str("");
  is >> nrewrt;
  ss << nrewrt;
  if (!(ss >> s.nrewrt)) {
    printIO("REPLICA", "NREWRT", nrewrt, "> 0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of READTRAJ block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, irottrans &s) {
  string e, rtc, rtclast;
  stringstream ss;
  s.found = 1;
  is >> rtc;
  ss << rtc;
  if (!(ss >> s.rtc)) {
    printIO("ROTTRANS", "RTC", rtc, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> rtclast;
  ss << rtclast;
  if (!(ss >> s.rtclast)) {
    printIO("ROTTRANS", "RTCLAST", rtclast, ">0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of ROTTRANS block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, istep &s) {
  string e, nstlim, t, dt;
  stringstream ss;
  s.found = 1;
  is >> nstlim;
  ss << nstlim;
  if (!(ss >> s.nstlim)) {
    printIO("STEP", "NSTLIM", nstlim, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> t;
  ss << t;
  if (!(ss >> s.t)) {
    printIO("STEP", "T", t, ">=0.0");
  }
  ss.clear();
  ss.str("");
  is >> dt;
  ss << dt;
  if (!(ss >> s.dt)) {
    printIO("STEP", "DT", dt, ">0.0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of STEP block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, istochdyn &s) {
  string e, ntsd, ntfr, nsfr, nbref, rcutf, cfric, tempsd;
  stringstream ss;
  s.found = 1;
  is >> ntsd;
  ss << ntsd;
  if (!(ss >> s.ntsd)) {
    printIO("STOCHDYN", "NTSD", ntsd, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> ntfr;
  ss << ntfr;
  if (!(ss >> s.ntfr)) {
    printIO("STOCHDYN", "NTFR", ntfr, "0..3");
  }
  ss.clear();
  ss.str("");
  is >> nsfr;
  ss << nsfr;
  if (!(ss >> s.nsfr)) {
    printIO("STOCHDYN", "NSFR", nsfr, ">0");
  }
  ss.clear();
  ss.str("");
  is >> nbref;
  ss << nbref;
  if (!(ss >> s.nbref)) {
    printIO("STOCHDYN", "NBREF", nbref, ">0");
  }
  ss.clear();
  ss.str("");
  is >> rcutf;
  ss << rcutf;
  if (!(ss >> s.rcutf)) {
    printIO("STOCHDYN", "RCUTF", rcutf, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> cfric;
  ss << cfric;
  if (!(ss >> s.cfric)) {
    printIO("STOCHDYN", "CFRIC", cfric, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> tempsd;
  ss << tempsd;
  if (!(ss >> s.tempsd)) {
    printIO("STOCHDYN", "TEMPSD", tempsd, ">=0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of STOCHDYN block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, isystem &s) {
  string e, npm, nsm;
  stringstream ss;
  s.found = 1;
  is >> npm;
  ss << npm;
  if (!(ss >> s.npm)) {
    printIO("SYSTEM", "NPM", npm, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> nsm;
  ss << nsm;
  if (!(ss >> s.nsm)) {
    printIO("SYSTEM", "NSM", nsm, ">=0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of SYSTEM block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ithermostat &s) {
  string e, ntt, s_ntbth, s_ntset, s_I, ntbtyp, tembth, s_ntbvar, s_taubth, ntscpl,
          ntstyp, ntscns, s_ntsgt, s_ntsgtg;
  stringstream ss;
  int ntbth, ntset, I, ntbvar, taubth, ntsgt, ntsgtg;
  s.found = 1;
  is >> ntt;
  ss << ntt;
  if (!(ss >> s.ntt)) {
    printIO("THERMOSTAT", "NTT", ntt, "0, 1");
  }
  ss.clear();
  ss.str("");
  is >> s_ntbth;
  ss << s_ntbth;
  if (!(ss >> ntbth)) {
    printIO("THERMOSTAT", "NTBTH", s_ntbth, ">= 0");
  }
  ss.clear();
  ss.str("");
  is >> s_ntset;
  ss << s_ntset;
  if (!(ss >> ntset)) {
    printIO("THERMOSTAT", "NTSET", s_ntset, ">= 0");
  }
  ss.clear();
  ss.str("");
  for (int i = 0; i < ntbth; ++i) {
    class ithermostat::tbath b;
    is >> s_I;
    ss << s_I;
    if(!(ss >> I)){
      stringstream msg;
      msg << "I(" << i+1 << ")";
      printIO("THERMOSTAT", msg.str(), s_I, "1..NTBATH");
    }
    ss.clear();
    ss.str("");
    is >> ntbtyp;
    ss << ntbtyp;
    if(!(ss >> b.ntbtyp)){
      stringstream msg;
      msg << "NTBTYP(" << i+1 << ")";
      printIO("THERMOSTAT", msg.str(), ntbtyp, "0..3");
    }
    ss.clear();
    ss.str("");
    is >> tembth;
    ss << tembth;
    if(!(ss >> b.tembth)){
      std::stringstream msg;
      msg << "TEMBTH(" << i+1 << ")";
      printIO("THERMOSTAT", msg.str(), ss.str(), ">= 0");
    }
    ss.clear();
    ss.str("");
    is >> s_ntbvar;
    ss << s_ntbvar;
    if(!(ss >> ntbvar)){
      std::stringstream msg;
      msg << "NTBVAR(" << i+1 << ")";
      printIO("THERMOSTAT", msg.str(), ss.str(), ">= 0 and <= MAX_NTBVAR");
    }
    ss.clear();
    ss.str("");
    for(int j = 0; j < ntbvar; ++j){
      is >> s_taubth;
      ss << s_taubth;
      if (!(ss >> taubth)){
        std::stringstream msg;
        msg << "TAUBTH(" << i+1 << ", " << j + 1 << ")";
        printIO("THERMOSTAT", msg.str(), ss.str(), ">= 0");
      }
      ss.clear();
      ss.str("");
      b.taubth.push_back(taubth);
    }
    s.baths.push_back(b);
  }
  for(int i = 0; i < ntset; ++i){
    class ithermostat::dofgroup d;
    is >> ntscpl;
    ss << ntscpl;
    if (!(ss >> d.ntscpl)){
      stringstream msg;
      msg << "NTSCPL(" << i+1 << ")";
      printIO("THERMOSTAT", msg.str(), ntscpl, "1..NTBTH");
    }
    ss.clear();;
    ss.str("");
    is >> ntstyp;
    ss << ntstyp;
    if(!(ss >> d.ntstyp)){
      stringstream msg;
      msg << "NTSTYP(" << i+1 << ")";
      printIO("THERMOSTAT", msg.str(), ntstyp, "0..2");
    }
    ss.clear();
    ss.str("");
    is >> ntscns;
    ss << ntscns;
    if(!(ss >> d.ntscns)){
      stringstream msg;
      msg << "NTSCNS(" << i+1 << ")";
      printIO("THERMOSTAT", msg.str(), ntscns, "0..2");
    }
    ss.clear();
    ss.str("");
    is >> s_ntsgt;
    ss << s_ntsgt;
    if(!(ss >> ntsgt)){
      stringstream msg;
      msg << "NTSGT(" << i+1 << ")";
      printIO("THERMOSTAT", msg.str(), s_ntsgt, ">= -2 and <=MAX_NTSGT");
    }
    if(ntsgt == -2){
      is >> s_ntsgtg;
      ss << s_ntsgtg;
      if (!(ss >> ntsgtg)) {
        stringstream msg;
        msg << "NTSGTG(" << i + 1 << ")";
        printIO("THERMOSTAT", msg.str(), s_ntsgt, ": see manual");
      }
      d.ntsgtg.push_back(ntsgtg);
      ss.clear();
      ss.str("");
      is >> s_ntsgtg;
      ss << s_ntsgtg;
      if (!(ss >> ntsgtg)) {
        stringstream msg;
        msg << "NTSGTG(" << i + 1 << ")";
        printIO("THERMOSTAT", msg.str(), s_ntsgt, ": see manual");
      }
      d.ntsgtg.push_back(ntsgtg);
      ss.clear();
      ss.str("");
    } else if (ntsgt > 1) {
      for (int j = 0; j < ntsgt; ++j) {
        is >> s_ntsgtg;
        ss << s_ntsgtg;
        if (!(ss >> ntsgtg)) {
          stringstream msg;
          msg << "NTSGTG(" << i + 1 << ", " << j + 1 << ")";
          printIO("THERMOSTAT", msg.str(), s_ntsgt, ": see manual");
        }
        d.ntsgtg.push_back(ntsgtg);
        ss.clear();
        ss.str("");
      }
    }
    s.dofgroups.push_back(d);
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of THERMOSTAT block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iumbrella &s) {
  string e, ntus, uscst1, uscst2, usref1, usref2;
  stringstream ss;
  s.found = 1;
  is >> ntus;
  ss << ntus;
  if (!(ss >> s.ntus)) {
    printIO("UMBRELLA", "NTUS", ntus, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> uscst1;
  ss << uscst1;
  if (!(ss >> s.uscst1)) {
    printIO("UMBRELLA", "USCST1", uscst1, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> uscst2;
  ss << uscst2;
  if (!(ss >> s.uscst2)) {
    printIO("UMBRELLA", "USCST2", uscst2, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> usref1;
  ss << usref1;
  if (!(ss >> s.usref1)) {
    printIO("UMBRELLA", "USREF1", usref1, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> usref2;
  ss << usref2;
  if (!(ss >> s.usref2)) {
    printIO("UMBRELLA", "USREF2", usref2, ">=0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of UMBRELLA block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ivirial &s) {
  string e, ntv, ntvg;
  stringstream ss;
  s.found = 1;
  is >> ntv;
  ss << ntv;
  if (!(ss >> s.ntv)) {
    printIO("VIRIAL", "NTV", ntv, "0,1");
  }
  ss.clear();
  ss.str("");
  is >> ntvg;
  ss << ntvg;
  if (!(ss >> s.ntvg)) {
    printIO("VIRIAL", "NTVG", ntvg, "0..3");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of VIRIAL block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, iwritetraj &s) {
  string e, ntwx, ntwse, ntwv, ntwf, ntwe, ntwg, ntwb;
  stringstream ss;
  s.found = 1;
  is >> ntwx;
  ss << ntwx;
  if (!(ss >> s.ntwx)) {
    printIO("WRITETRAJ", "NTWX", ntwx, "an integer.");
  }
  ss.clear();
  ss.str("");
  is >> ntwse;
  ss << ntwse;
  if (!(ss >> s.ntwse)) {
    printIO("WRITETRAJ", "NTWSE", ntwse, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> ntwv;
  ss << ntwv;
  if (!(ss >> s.ntwv)) {
    printIO("WRITETRAJ", "NTWV", ntwv, "an integer.");
  }
  ss.clear();
  ss.str("");
  is >> ntwf;
  ss << ntwf;
  if (!(ss >> s.ntwf)) {
    printIO("WRITETRAJ", "NTWF", ntwf, "an integer.");
  }
  ss.clear();
  ss.str("");
  is >> ntwe;
  ss << ntwe;
  if (!(ss >> s.ntwe)) {
    printIO("WRITETRAJ", "NTWE", ntwe, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> ntwg;
  ss << ntwg;
  if (!(ss >> s.ntwg)) {
    printIO("WRITETRAJ", "NTWG", ntwg, ">=0");
  }
  ss.clear();
  ss.str("");
  is >> ntwb;
  ss << ntwb;
  if (!(ss >> s.ntwb)) {
    printIO("WRITETRAJ", "NTWB", ntwb, ">=0");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of WRITETRAJ block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ixrayres &s) {
  string e, ntxr, ntxle, cxr, ntwxr, ntwde, mtwde, ntwxm, cxtau, rdavg;
  stringstream ss;
  s.found = 1;
  is >> ntxr;
  ss << ntxr;
  if (!(ss >> s.ntxr)) {
    printIO("XRAYRES", "NTXR", ntxr, "NTXR must be -2..3");
  }
  ss.clear();
  ss.str("");
  is >> ntxle;
  ss << ntxle;
  if (!(ss >> s.ntxle)) {
    printIO("XRAYRES", "NTXLE", ntxle, "NTXLE must be 0,1");
  }
  ss.clear();
  ss.str("");
  is >> cxr;
  ss << cxr;
  if (!(ss >> s.cxr)) {
    printIO("XRAYRES", "CXR", cxr, "CXR must be >= 0.0");
  }
  ss.clear();
  ss.str();
  is >> ntwxr;
  ss << ntwxr;
  if (!(ss >> s.ntwxr)) {
    printIO("XRAYRES", "NTWXR", ntwxr, "NTWXR must be >= 0");
  }
  ss.clear();
  ss.str("");
  is >> ntwde;
  ss << ntwde;
  if (!(ss >> s.ntwde)) {
    printIO("XRAYRES", "NTWDE", ntwde, "NTWDE must be 0..3");
  }
  ss.clear();
  ss.str("");
  is >> ntwxm;
  ss << ntwxm;
  if (!(ss >> s.ntwxm)) {
    printIO("XRAYRES", "NTWXM", ntwxm, "NTWXM must be >= 0");
  }
  ss.clear();
  ss.str("");
  is >> cxtau;
  ss << cxtau;
  if (!(ss >> s.cxtau)) {
    printIO("XRAYRES", "CXTAU", cxtau, "CXTAU must be > 0.0");
  }
  ss.clear();
  ss.str("");
  is >> rdavg;
  ss << rdavg;
  if (!(ss >> s.rdavg)) {
    printIO("XRAYRES", "RDAVG", rdavg, "RDAVG must be 0,1");
  }
  is >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of XRAYRES block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
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
    if (buffer.size()) {


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
          cout << "Don't know anything about block " << buffer[0]
                  << ". Just storing data.\n";
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

    if (first == "GENBOX") {
      s.box.boxformat() = gcore::Box::genbox;
      int ntb;
      gmath::Vec k, l, m;
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
      s.box.K() = k;
      s.box.L() = l;
      s.box.M() = m;

    }

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
            << "#      NTT     NTBTH     NTSET\n"
            << setw(10) << gin.thermostat.ntt
            << setw(10) << gin.thermostat.baths.size()
            << setw(10) << gin.thermostat.dofgroups.size()
            << "\n"
            << "# I = 1 ... NTBTH\n"
            << "# I TEMBTH(I)  NTBVAR(I)   TAUBTH(I,1...NTVAR)\n";
    for (unsigned int i = 0; i < gin.thermostat.baths.size(); ++i) {
      os << setw(3) << i + 1
              << setw(10) << gin.thermostat.baths[i].tembth
              << setw(10) << gin.thermostat.baths[i].taubth.size();

      for (unsigned int j = 0; j < gin.thermostat.baths[i].taubth.size(); ++j) {
        os << setw(7) << gin.thermostat.baths[i].taubth[j];
      }
      os << "\n";
    }
    os << "# J = 1 ... NTGRP\n"
            << "# NTSCPL(J) NTSTYP(J) NTSCNS(J)  NTSGT(J)  NTSGTG(J,1...NTGT(J))\n";
    for (unsigned int j = 0; j < gin.thermostat.dofgroups.size(); ++j) {
      os << setw(10) << gin.thermostat.dofgroups[j].ntscpl
              << setw(10) << gin.thermostat.dofgroups[j].ntstyp
              << setw(10) << gin.thermostat.dofgroups[j].ntscns
              << setw(10) << gin.thermostat.dofgroups[j].ntsgt;
      if (gin.thermostat.dofgroups[j].ntsgt > 0) {
        for (unsigned int k = 0;
                k < gin.thermostat.dofgroups[j].ntsgtg.size(); ++k) {
          os << setw(7) << gin.thermostat.dofgroups[j].ntsgtg[k];
        }
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
    if (gin.constraint.ntcs == 1 || gin.constraint.ntcs == 2) {
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
            << "#     COS    EFIELD    MINFIELD  DAMP     WRITE\n"
            << setw(10) << gin.polarise.cos
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
            << "\nEND\n";
  }

  // PAIRLIST GENERATION

  // NEIGHBOURLIST (promd)
  if (gin.neighbourlist.found) {
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
            << "#    NTLES    NLEPOT    NTLESA\n"
            << setw(10) << gin.localelev.ntles
            << setw(10) << gin.localelev.nlepot
            << setw(10) << gin.localelev.ntlesa << endl
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
            << "# NTPW = 0 : binary\n"
            << "# NTPW = 1 : formatted\n"
            << "# NTWSE = configuration selection parameter\n"
            << "# =0: write normal trajectory\n"
            << "# >0: chose min energy for writing configurations\n"
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
