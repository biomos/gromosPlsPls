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
    for(int i = 0; i < 0; ++i) {
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
  int found, ntxr, ntwxr, ntwde, ntwxm, rdavg;
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
  if(!(ss >> s.ntp)) {
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
  if(e != "") {
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
  if(!(ss >> s.ntb)){
    printIO("BOUNDCOND", "NTB", ntb, "-1..2");
  }
  ss.clear();
  is >> ndfmin;
  ss << ndfmin;
  if(!(ss >> s.ndfmin)) {
    printIO("BOUNDCOND", "NDFMIN", ndfmin, ">= 0");
  }
  is >> e;
  if(e != "") {
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
  if(!(ss >> s.ntcgran)) {
    printIO("CGRAIN", "NTCGRAN", ntcgran, "0..2");
  }
  ss.clear();
  is >> eps;
  ss << eps;
  if(!(ss >> s.eps)){
    printIO("CGRAIN", "EPS", eps, ">= 0.0");
  }
  is >> e;
  if(e != "") {
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
  if(!(ss >> s.nscm)){
    stringstream msg;
    msg << "Could not convert \"" << nscm << "\" to an integer value in COMTRANSROT block";
    printError(msg.str());
  }
  is >> e;
  if(e != "") {
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
  if(!(ss >> s.ntchk)){
    printIO("CONSISTENCYCHECK", "NTCHK", ntchk, "0, 1");
  }
  ss.clear();
  is >> ntckf;
  ss << ntckf;
  if (!(ss >> s.ntckf)){
    printIO("CONSISTENCYCHECK", "NTCKF", ntckf, "0, 1");
  }
  ss.clear();
  is >> fdckf;
  ss << fdckf;
  if(!(ss >> s.fdckf)){
    printIO("CONSISTENCYCHECK", "FDCKF", fdckf, "> 0.0");
  }
  ss.clear();
  is >> ntckv;
  ss << ntckv;
  if(!(ss >> s.ntckv)){
    printIO("CONSISTENCYCHECK", "NTCKV", ntckv, "> 0.0");
  }
  ss.clear();
  is >> fdckv;
  ss << fdckv;
  if(!(ss >> s.fdckv)){
    printIO("CONSISTENCYCHECK", "FDCKV", fdckv, "> 0.0");
  }
  ss.clear();
  is >> ntckt;
  ss << ntckt;
  if (!(ss >> s.ntckt)){
    printIO("CONSISTENCYCHECK", "NTCKT", ntckt, "> 0.0");
  }
  ss.clear();
  is >> ntcke;
  ss << ntcke;
  if(!(ss >> s.ntcke)){
    printIO("CONSISTENCYCHECK", "NTCKE", ntcke, "> 0.0");
  }
  ss.clear();
  is >> ntckr;
  ss << ntckr;
  if(!(ss >> s.ntckr)){
    printIO("CONSISTENCYCHECK", "NTCKR", ntckr, "> 0.0");
  }
  ss.clear();
  is >> ntckl;
  ss << ntckl;
  if(!(ss >> s.ntckl)){
    printIO("CONSISTENCYCHECK", "NTCKL", ntckl, "> 0.0");
  }
  ss.clear();
  is >> fdckl;
  ss << fdckl;
  if(!(ss >> s.fdckl)){
    printIO("CONSISTENCYCHECK", "FDCKL", fdckl, "> 0.0");
  }
  ss.clear();
  is >> s_nackf;
  ss << s_nackf;
  if(!(ss >> nackf)) {
    printIO("CONSISTENCYCHECK", "NACKF", s_nackf, ">= 0");
  }
  ss.clear();
  if (nackf < 0) {
    std::stringstream ss;
    ss << nackf;
    printIO("CONSISTENCYCHECK", "NACKF", ss.str(), ">= 0");
  }
  for (int i = 0; i < nackf; i++) {
    is >> s_nckf;
    ss << s_nckf;
    if(!(ss >> nckf)){
      printIO("CONSISTENCYCHECK", "NCKF", s_nckf, ">= 1");
    }
    ss.clear();
    s.nckf.push_back(nckf);
  }
  is >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of CONSISTENCYCHECK block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iconstraint &s) {
  string e;
  s.found = 1;
  is >> s.ntc >> s.ntcp >> s.ntcp0[0];
  if (s.ntcp == 3) {
    is >> s.ntcp0[1] >> s.ntcp0[2];
  }
  is >> s.ntcs;
  if (s.ntcs != 4) {
    is >> s.ntcs0[0];
  }
  if (s.ntcs == 3) {
    is >> s.ntcs0[1] >> s.ntcs0[2];
  }
  is >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of CONSTRAINT block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, icovalentform &s) {
  string e;
  s.found = 1;
  is >> s.ntbbh >> s.ntbah >> s.ntbdn >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of COVALENTFORM block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, idebug &s) {
  string e;
  int nrd;
  s.found = 1;
  is >> nrd;
  if (nrd < 0) {
    std::stringstream ss;
    ss << nrd;
    printIO("DEBUG", "NRD", ss.str(), ">= 0");
  }
  for (int i = 0; i < nrd; ++i) {
    class idebug::routine r;
    is >> r.piider;
    is >> r.iiideo;
    s.routines.push_back(r);
  }
  is >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of DEBUG block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, idihedralres &s) {
  string e;
  s.found = 1;
  is >> s.ntdlr >> s.cdlr >> s.philin >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of DIHEDRALRES block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, idistanceres &s) {
  string e;
  s.found = 1;
  is >> s.ntdir >> s.ntdira >> s.cdir >> s.dir0 >> s.taudir >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of DISTANCERES block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ienergymin &s) {
  string e;
  s.found = 1;
  is >> s.ntem >> s.ncyc >> s.dele >> s.dx0 >> s.dxm >> s.nmin >> s.flim >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of ENERGYMIN block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, iewarn &s) {
  string e;
  s.found = 1;
  is >> s.maxener >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of EWARN block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, iforce &s) {
  string e;
  s.found = 1;
  int negr, nre;
  for (int i = 0; i < 10; i++) {
    is >> s.ntf[i];
  }
  is >> negr;
  if (negr <= 0) {
    std::stringstream ss;
    ss << negr;
    printIO("FORCE", "NEGR", ss.str(), "> 0");
  }
  for (int i = 0; i < negr; i++) {
    is >> nre;
    s.nre.push_back(nre);
  }
  is >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of FORCE block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, igeomconstraints &s) {
  string e;
  s.found = 1;
  is >> s.ntcph >> s.ntcpn >> s.ntcs >> s.shktol >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of GEOMCONSTRAINTS block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, igromos96compat &s) {
  string e;
  s.found = 1;
  is >> s.ntnb96 >> s.ntr96 >> s.ntp96 >> s.ntg96 >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of GROMOS96COMPAT block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, iinitialise &s) {
  string e;
  s.found = 1;
  is >> s.ntivel >> s.ntishk >> s.ntinht >> s.ntinhb;
  is >> s.ntishi >> s.ntirtc >> s.nticom;
  is >> s.ntisti >> s.ig >> s.tempi >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of INITIALISE block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, iinnerloop &s) {
  string e;
  s.found = 1;
  is >> s.ntilm >> s.ntils >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of INNERLOOP block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iintegrate &s) {
  string e;
  s.found = 1;
  is >> s.nint >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of INTEGRATE block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ijvalueres &s) {
  string e;
  s.found = 1;
  is >> s.ntjvr >> s.ntjvra >> s.cjvr >> s.taujvr >> s.le
          >> s.ngrid >> s.delta >> s.write >> e;
  if(e != "") {
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
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of LAMBDAS block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ilocalelev &s) {
  string e;
  s.found = 1;
  is >> s.ntles >> s.nlepot >> s.ntlesa;
  if (s.nlepot < 0) {
    std::stringstream ss;
    ss << s.nlepot;
    printIO("LOCALELEV", "NLEPOT", ss.str(), ">= 0");
  }
  int id, erfl;
  for (int i = 0; i < s.nlepot; i++) {
    is >> id >> erfl;
    s.nlepid_ntlerf.insert(pair<int, int> (id, erfl));
  }
  if (int(s.nlepid_ntlerf.size()) != s.nlepot) {
    std::stringstream ss, si;
    ss << s.nlepid_ntlerf.size() << " potentials";
    si << s.nlepot;
    printIO("LOCALELEV", "NLEPID and NTLEFR", ss.str(), si.str());
  }
  is >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of LOCALELEV block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, imultibath &s) {
  string e;
  double temp0, tau;
  int last, combath, irbath;
  s.found = 1;
  is >> s.algorithm;
  is >> s.nbaths;
  if (s.nbaths < 0) {
    std::stringstream ss;
    ss << s.nbaths;
    printIO("MULTIBATH", "NBATHS", ss.str(), "> 0");
  }
  for (int i = 0; i < s.nbaths; ++i) {
    is >> temp0 >> tau;
    s.temp0.push_back(temp0);
    s.tau.push_back(tau);
  }
  is >> s.dofset;
  if (s.dofset < 0) {
    std::stringstream ss;
    ss << s.dofset;
    printIO("MULTIBATH", "DOFSET", ss.str(), ">= 0");
  }
  for (int i = 0; i < s.dofset; ++i) {
    is >> last >> combath >> irbath;
    s.last.push_back(last);
    s.combath.push_back(combath);
    s.irbath.push_back(irbath);
  }
  is >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of MULTIBATH block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, imulticell &s) {
  string e;
  s.found = 1;
  is >> s.ntm >> s.ncella >> s.ncellb >> s.ncellc
          >> s.tolpx >> s.tolpv >> s.tolpf >> s.tolpfw >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of MULTICELL block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, imultistep &s) {
  string e;
  s.found = 1;
  is >> s.steps >> s.boost >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of MULTISTEP block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ineighbourlist &s) {
  string e;
  s.found = 1;
  is >> s.nmprpl >> s.nuprpl >> s.rcprpl >> s.grprpl >> s.nmtwpl;
  is >> s.nutwpl >> s.rstwpl >> s.rltwpl >> s.nuirin >> s.nusrin;
  is >> s.nmtwin >> s.rctwin >> s.ncgcen >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of NEIGBOURLIST block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, inonbonded &s) {
  string e;
  s.found = 1;
  is >> s.nlrele;
  is >> s.appak >> s.rcrf >> s.epsrf;
  is >> s.nshape >> s.ashape >> s.na2clc >> s.tola2 >> s.epsls;
  is >> s.nkx >> s.nky >> s.nkz >> s.kcut;
  is >> s.ngx >> s.ngy >> s.ngz >> s.nasord >> s.nfdord >> s.nalias
          >> s.nspord;
  is >> s.nqeval >> s.faccur >> s.nrdgrd >> s.nwrgrd;
  is >> s.nlrlj >> s.slvdns >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of NONBONDED block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ioveralltransrot &s) {
  string e;
  s.found = 1;
  is >> s.ncmtr >> s.ncmro >> s.cmamx >> s.cmamy >> s.cmamz >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of OVERALLTRANSROT block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ipairlist &s) {
  string e, alg, siz, typ;
  s.found = 1;
  is >> s.algorithm >> s.nsnb >> s.rcutp >> s.rcutl >> s.size >> s.type >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of PAIRLIST block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ipathint &s) {
  string e;
  s.found = 1;
  is >> s.ntpi >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of PATHINT block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iperscale &s) {
  string e;
  s.found = 1;
  is >> s.restype >> s.kdih >> s.kj >> s.t >> s.diff >> s.ratio >> s.read >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of PERSCALE block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iperturbation &s) {
  string e;
  s.found = 1;
  is >> s.ntg >> s.nrdgl >> s.rlam >> s.dlamt
          >> s.alphlj >> s.alphc >> s.nlam >> s.nscale >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of PERTURBATION block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ipolarise &s) {
  string e;
  s.found = 1;
  is >> s.cos >> s.efield >> s.minfield >> s.damp >> s.write >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of POLARISE block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ipositionres &s) {
  string e;
  s.found = 1;
  is >> s.ntpor >> s.ntporb >> s.ntpors >> s.cpor >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of POSITIONRES block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ipressurescale &s) {
  string e;
  s.found = 1;
  is >> s.couple >> s.scale >> s.comp >> s.taup >> s.virial;
  is >> s.pres0[0][0] >> s.pres0[0][1] >> s.pres0[0][2]
	  >> s.pres0[1][0] >> s.pres0[1][1] >> s.pres0[1][2]
	  >> s.pres0[2][0] >> s.pres0[2][1] >> s.pres0[2][2] >> e;
  if (e != "") {
    stringstream ss;
    ss << "unexpected end of PRESSURESCALE block, read \"" << e << "\" instead of \"END\"";
    printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iprintout &s) {
  string e;
  s.found = 1;
  is >> s.ntpr >> s.ntpp >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of PRINTOUT block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, irandomnumbers &s) {
  string e;
  s.found = 1;
  is >> s.ntrng >> s.ntgsl >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of RANDOMNUMBERS block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, ireadtraj &s) {
  string e;
  s.found = 1;
  is >> s.ntrd >> s.ntrn >> s.ntrb >> s.ntshk >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of READTRAJ block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ireplica &s) {
  string e;
  int nret, nrelam;
  double dum;
  s.found = 1;
  is >> nret;
  if (nret < 0) {
    std::stringstream ss;
    ss << nret;
    printIO("REPLICA", "NRET", ss.str(), ">= 0");
  }
  for (int i = 0; i < nret; ++i) {
    is >> dum;
    s.ret.push_back(dum);
  }
  is >> s.lrescale;
  is >> nrelam;
  if (nrelam < 0) {
    std::stringstream ss;
    ss << nrelam;
    printIO("REPLICA", "NRELAM", ss.str(), ">= 0");
  }
  for (int i = 0; i < nrelam; ++i) {
    is >> dum;
    s.relam.push_back(dum);
    is >> dum;
    s.rets.push_back(dum);
  }
  is >> s.nretrial >> s.nrequil >> s.nrejob >> s.nrewrt >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of READTRAJ block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, irottrans &s) {
  string e;
  s.found = 1;
  is >> s.rtc >> s.rtclast >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of ROTTRANS block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, istep &s) {
  string e;
  s.found = 1;
  is >> s.nstlim >> s.t >> s.dt >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of STEP block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, istochdyn &s) {
  string e;
  s.found = 1;
  is >> s.ntsd >> s.ntfr >> s.nsfr >> s.nbref >> s.rcutf >> s.cfric >> s.tempsd >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of STOCHDYN block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, isystem &s) {
  string e;
  s.found = 1;
  is >> s.npm >> s.nsm >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of SYSTEM block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ithermostat &s) {
  string e;
  int idum, dum, ntbth, ntset;
  double tau;
  s.found = 1;
  is >> s.ntt >> ntbth >> ntset;
  if (ntbth < 0) {
    std::stringstream ss;
    ss << ntbth;
    printIO("THERMOSTAT", "NTBTH", ss.str(), ">=0");
  }
  if (ntset < 0) {
    std::stringstream ss;
    ss << ntset;
    printIO("THERMOSTAT", "NTSET", ss.str(), ">= 0");
  }
  for (int i = 0; i < ntbth; ++i) {
    class ithermostat::tbath b;
    is >> idum >> b.ntbtyp >> b.tembth >> dum;
    if (idum != i + 1) {
      std::stringstream ss, si, sj;
      ss << idum;
      si << "I[" << i + 1 << "]";
      sj << i + 1;
      printIO("THERMOSTAT", si.str(), ss.str(), sj.str());
    }
    if (dum < 0) {
      std::stringstream ss, si;
      ss << dum;
      si << "NTBVAR[" << i + 1 << "]";
      printIO("THERMOSTAT", si.str(), ss.str(), ">= 0.0");
    }
    for (int j = 0; j < dum; ++j) {
      is >> tau;
      b.taubth.push_back(tau);
    }
    s.baths.push_back(b);
  }
  for (int i = 0; i < ntset; ++i) {
    class ithermostat::dofgroup d;
    is >> d.ntscpl >> d.ntstyp >> d.ntscns >> d.ntsgt;
    int gt = 0, ntsgtg;
    if (d.ntsgt <= 0) gt = 0;
    else gt = d.ntsgt;
    for (int k = 0; k < gt; k++) {
      is >> ntsgtg;
      d.ntsgtg.push_back(ntsgtg);
    }
    s.dofgroups.push_back(d);
  }
  is >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of THERMOSTAT block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }
  return is;
}

istringstream & operator>>(istringstream &is, iumbrella &s) {
  string e;
  s.found = 1;
  is >> s.ntus >> s.uscst1 >> s.uscst2 >> s.usref1 >> s.usref2 >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of UMBRELLA block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ivirial &s) {
  string e;
  s.found = 1;
  is >> s.ntv >> s.ntvg >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of VIRIAL block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, iwritetraj &s) {
  string e;
  s.found = 1;
  is >> s.ntwx >> s.ntwse >> s.ntwv >> s.ntwf >> s.ntwe >> s.ntwg >> s.ntwb >> e;
  if(e != "") {
      stringstream ss;
      ss << "unexpected end of WRITETRAJ block, read \"" << e << "\" instead of \"END\"";
      printError(ss.str());
  }

  return is;
}

istringstream & operator>>(istringstream &is, ixrayres &s) {
  string e;
  s.found = 1;
  is >> s.ntxr >> s.cxr >> s.ntwxr >> s.ntwde >> s.ntwxm >> s.cxtau >> s.rdavg >> e;
  if(e != "") {
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
              << setw(10) << gin.multibath.tau[i];
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
            << "#    NTXR   CXR   NTWXR   NTWDE   NTWXM   CXTAU  RDAVG\n"
            << setw(10) << gin.xrayres.ntxr
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
